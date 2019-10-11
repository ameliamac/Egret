#  ___________________________________________________________________________
#
#  EGRET: Electrical Grid Research and Engineering Tools
#  Copyright 2019 National Technology & Engineering Solutions of Sandia, LLC
#  (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S.
#  Government retains certain rights in this software.
#  This software is distributed under the Revised BSD License.
#  ___________________________________________________________________________

"""
This module provides functions that create the modules an ACOPF SOCP relaxation

"""
import pyomo.environ as pe
import egret.model_library.transmission.tx_utils as tx_utils
import egret.model_library.transmission.tx_calc as tx_calc
import operator as op
import egret.model_library.transmission.bus as libbus
import egret.model_library.transmission.branch as libbranch
import egret.model_library.transmission.gen as libgen
import math 

from egret.model_library.defn import FlowType, CoordinateType
from egret.data.model_data import map_items, zip_items


    
def create_socp_acopf_model(model_data, include_feasibility_slack=False):
    md = model_data.clone_in_service()
    tx_utils.scale_ModelData_to_pu(md, inplace = True)

    gens = dict(md.elements(element_type='generator'))
    buses = dict(md.elements(element_type='bus'))
    branches = dict(md.elements(element_type='branch'))
    loads = dict(md.elements(element_type='load'))
    shunts = dict(md.elements(element_type='shunt'))

    gen_attrs = md.attributes(element_type='generator')
    bus_attrs = md.attributes(element_type='bus')
    branch_attrs = md.attributes(element_type='branch')
    load_attrs = md.attributes(element_type='load')
    shunt_attrs = md.attributes(element_type='shunt')

    inlet_branches_by_bus, outlet_branches_by_bus = \
        tx_utils.inlet_outlet_branches_by_bus(branches, buses)
    gens_by_bus = tx_utils.gens_by_bus(buses, gens)

    model = pe.ConcreteModel()

    ### declare (and fix) the loads at the buses
    bus_p_loads, bus_q_loads = tx_utils.dict_of_bus_loads(buses, loads)

    libbus.declare_var_pl(model, bus_attrs['names'], initialize=bus_p_loads)
    libbus.declare_var_ql(model, bus_attrs['names'], initialize=bus_q_loads)
    model.pl.fix()
    model.ql.fix()

    ### declare the fixed shunts at the buses
    bus_bs_fixed_shunts, bus_gs_fixed_shunts = tx_utils.dict_of_bus_fixed_shunts(buses, shunts)


    ### declare the rectangular voltages
    neg_v_max = map_items(op.neg, bus_attrs['v_max'])
    vr_init = {k: bus_attrs['vm'][k] * pe.cos(bus_attrs['va'][k]) for k in bus_attrs['vm']}
    # libbus.declare_var_vr(model, bus_attrs['names'], initialize=vr_init,
    #                       bounds=zip_items(neg_v_max, bus_attrs['v_max'])
    #                       )

    vj_init = {k: bus_attrs['vm'][k] * pe.sin(bus_attrs['va'][k]) for k in bus_attrs['vm']}
    # libbus.declare_var_vj(model, bus_attrs['names'], initialize=vj_init,
    #                       bounds=zip_items(neg_v_max, bus_attrs['v_max'])
    #                       )


    # w variable for socp
    vj2_init = {k: (bus_attrs['vm'][k] * pe.sin(bus_attrs['va'][k]))**2 for k in bus_attrs['vm']}
    vr2_init = {k: (bus_attrs['vm'][k] * pe.cos(bus_attrs['va'][k]))**2 for k in bus_attrs['vm']}
    w_init = {k: vr2_init[k]+vj2_init[k] for k in vj2_init}
    wub = {k:bus_attrs['v_max'][k]**2 for k in bus_attrs['v_max']}
    wlb = {k:bus_attrs['v_min'][k]**2 for k in bus_attrs['v_min']}
    ## v_min**2 <= w <= v_max**2
    libbus.declare_var_w(model, bus_attrs['names'], initialize = w_init, 
                        bounds =zip_items(wlb,wub))



    ### include the feasibility slack for the bus balances
    p_rhs_kwargs = {}
    q_rhs_kwargs = {}
    if include_feasibility_slack:
        p_rhs_kwargs, q_rhs_kwargs, penalty_expr = _include_feasibility_slack(model, bus_attrs, gen_attrs, bus_p_loads, bus_q_loads)

    ### fix the reference bus
    # ref_bus = md.data['system']['reference_bus']
    # ref_angle = md.data['system']['reference_bus_angle']
    # if ref_angle != 0.0:
    #     libbus.declare_eq_ref_bus_nonzero(model, ref_angle, ref_bus)
    # else:
    #     model.vj[ref_bus].fix(0.0)
    #     model.vr[ref_bus].setlb(0.0)

    ### declare the generator real and reactive power
    pg_init = {k: (gen_attrs['p_min'][k] + gen_attrs['p_max'][k]) / 2.0 for k in gen_attrs['pg']}
    libgen.declare_var_pg(model, gen_attrs['names'], initialize=pg_init,
                          bounds=zip_items(gen_attrs['p_min'], gen_attrs['p_max'])
                          )

    qg_init = {k: (gen_attrs['q_min'][k] + gen_attrs['q_max'][k]) / 2.0 for k in gen_attrs['qg']}
    libgen.declare_var_qg(model, gen_attrs['names'], initialize=qg_init,
                          bounds=zip_items(gen_attrs['q_min'], gen_attrs['q_max'])
                          )

    ### declare the current flows in the branches
    s_max = {k: branches[k]['rating_long_term'] for k in branches.keys()}
    s_lbub = dict()
    for k in branches.keys():
        if s_max[k] is None:
            s_lbub[k] = (None, None)
        else:
            s_lbub[k] = (-s_max[k],s_max[k])
    pf_bounds = s_lbub
    pt_bounds = s_lbub
    qf_bounds = s_lbub
    qt_bounds = s_lbub
    pf_init = dict()
    pt_init = dict()
    qf_init = dict()
    qt_init = dict()
    cbk_init = dict()
    sbk_init = dict()
    cbk_bounds = dict()
    sbk_bounds = dict()
    for branch_name, branch in branches.items():
        from_bus = branch['from_bus']
        to_bus = branch['to_bus']
        ba_max = branch['angle_diff_max']
        ba_min = branch['angle_diff_min']
        y_matrix = tx_calc.calculate_y_matrix_from_branch(branch)
        ifr_init = tx_calc.calculate_ifr(vr_init[from_bus], vj_init[from_bus], vr_init[to_bus],
                                         vj_init[to_bus], y_matrix)
        ifj_init = tx_calc.calculate_ifj(vr_init[from_bus], vj_init[from_bus], vr_init[to_bus],
                                         vj_init[to_bus], y_matrix)
        itr_init = tx_calc.calculate_itr(vr_init[from_bus], vj_init[from_bus], vr_init[to_bus],
                                         vj_init[to_bus], y_matrix)
        itj_init = tx_calc.calculate_itj(vr_init[from_bus], vj_init[from_bus], vr_init[to_bus],
                                         vj_init[to_bus], y_matrix)
        pf_init[branch_name] = tx_calc.calculate_p(ifr_init, ifj_init, vr_init[from_bus], vj_init[from_bus])
        pt_init[branch_name] = tx_calc.calculate_p(itr_init, itj_init, vr_init[to_bus], vj_init[to_bus])
        qf_init[branch_name] = tx_calc.calculate_q(ifr_init, ifj_init, vr_init[from_bus], vj_init[from_bus])
        qt_init[branch_name] = tx_calc.calculate_q(itr_init, itj_init, vr_init[to_bus], vj_init[to_bus])

        #SOCP related variable bounds and init
        cbk_init[from_bus,to_bus] = vr_init[from_bus]*vr_init[to_bus] + vj_init[from_bus]*vj_init[to_bus]
        sbk_init[from_bus,to_bus] = vr_init[from_bus]*vj_init[to_bus] - vr_init[to_bus]*vj_init[from_bus]

        if ba_max is None and ba_min is None:
          cbk_bounds[from_bus,to_bus] = ( bus_attrs['v_min'][from_bus]*bus_attrs['v_min'][to_bus]*min(pe.cos(-math.pi/2),pe.cos(math.pi/2)),
                                      bus_attrs['v_max'][from_bus]*bus_attrs['v_max'][to_bus]*1.0)
          sbk_bounds[from_bus,to_bus] = ( bus_attrs['v_max'][from_bus]*bus_attrs['v_max'][to_bus]*pe.sin(-math.pi/2),
                                      bus_attrs['v_max'][from_bus]*bus_attrs['v_max'][to_bus]*pe.sin(math.pi/2))
        if ba_max > 0 and ba_min < 0:
          cbk_bounds[from_bus,to_bus] = (bus_attrs['v_min'][from_bus]*bus_attrs['v_min'][to_bus]*min(pe.cos(ba_max * math.pi/180),pe.cos(ba_min * math.pi/180)),
                                      bus_attrs['v_max'][from_bus]*bus_attrs['v_max'][to_bus]*1.0)
          sbk_bounds[from_bus,to_bus] = ( bus_attrs['v_max'][from_bus]*bus_attrs['v_max'][to_bus]*pe.sin(ba_min * math.pi/180),
                                      bus_attrs['v_max'][from_bus]*bus_attrs['v_max'][to_bus]*pe.sin(math.pi/2))
        if ba_max <= 0:
          cbk_bounds[from_bus,to_bus] = (bus_attrs['v_min'][from_bus]*bus_attrs['v_min'][to_bus]*pe.cos(ba_min * math.pi/180),
                                      bus_attrs['v_max'][from_bus]*bus_attrs['v_max'][to_bus]*pe.cos(ba_max * math.pi/180))
          sbk_bounds[from_bus,to_bus] = ( bus_attrs['v_max'][from_bus]*bus_attrs['v_max'][to_bus]*pe.sin(ba_min * math.pi/180),
                                      bus_attrs['v_max'][from_bus]*bus_attrs['v_max'][to_bus]*pe.sin(math.pi/2))
        if ba_min >= 0:
          cbk_bounds[from_bus,to_bus] = (bus_attrs['v_min'][from_bus]*bus_attrs['v_min'][to_bus]*pe.cos(ba_max * math.pi/180),
                                      bus_attrs['v_max'][from_bus]*bus_attrs['v_max'][to_bus]*pe.cos(ba_min * math.pi/180))
          sbk_bounds[from_bus,to_bus] = ( bus_attrs['v_max'][from_bus]*bus_attrs['v_max'][to_bus]*pe.sin(ba_min * math.pi/180),
                                      bus_attrs['v_max'][from_bus]*bus_attrs['v_max'][to_bus]*pe.sin(math.pi/2))

    
    #print(bus_attrs)
    #print(branch_attrs)

    libbranch.declare_var_pf(model=model,
                             index_set=branch_attrs['names'],
                             initialize=pf_init,
                             bounds=pf_bounds
                             )
    libbranch.declare_var_pt(model=model,
                             index_set=branch_attrs['names'],
                             initialize=pt_init,
                             bounds=pt_bounds
                             )
    libbranch.declare_var_qf(model=model,
                             index_set=branch_attrs['names'],
                             initialize=qf_init,
                             bounds=qf_bounds
                             )
    libbranch.declare_var_qt(model=model,
                             index_set=branch_attrs['names'],
                             initialize=qt_init,
                             bounds=qt_bounds
                             )


    bus_pairs = zip_items(branch_attrs['from_bus'],branch_attrs['to_bus'])
    unique_bus_pairs = list(set([val for idx,val in bus_pairs.items()]))

    libbranch.declare_var_c(model = model,
                             index_set = unique_bus_pairs,
                             initialize = cbk_init,
                             bounds = cbk_bounds
                             )

    libbranch.declare_var_s(model = model,
                             index_set = unique_bus_pairs,
                             initialize = sbk_init,
                             bounds = sbk_bounds
                             )

    ### declare the branch power flow constraints
    libbranch.declare_eq_branch_power_socp(model=model,
                                      index_set=branch_attrs['names'],
                                      branches=branches,
                                      branch_attrs=branch_attrs,
                                      coordinate_type=CoordinateType.RECTANGULAR
                                      )

    ### declare the pq balances
    libbus.declare_eq_p_balance_socp(model=model,
                                index_set=bus_attrs['names'],
                                bus_p_loads=bus_p_loads,
                                gens_by_bus=gens_by_bus,
                                bus_gs_fixed_shunts=bus_gs_fixed_shunts,
                                inlet_branches_by_bus=inlet_branches_by_bus,
                                outlet_branches_by_bus=outlet_branches_by_bus,
                                coordinate_type=CoordinateType.RECTANGULAR,
                                **p_rhs_kwargs
                                )

    libbus.declare_eq_q_balance_socp(model=model,
                                index_set=bus_attrs['names'],
                                bus_q_loads=bus_q_loads,
                                gens_by_bus=gens_by_bus,
                                bus_bs_fixed_shunts=bus_bs_fixed_shunts,
                                inlet_branches_by_bus=inlet_branches_by_bus,
                                outlet_branches_by_bus=outlet_branches_by_bus,
                                coordinate_type=CoordinateType.RECTANGULAR,
                                **q_rhs_kwargs
                                )

    ### declare the thermal limits
    libbranch.declare_ineq_s_branch_thermal_limit(model=model,
                                                  index_set=branch_attrs['names'],
                                                  branches=branches,
                                                  s_thermal_limits=s_max,
                                                  flow_type=FlowType.POWER
                                                  )

    ### declare the voltage min and max inequalities
    # libbus.declare_ineq_vm_bus_lbub(model=model,
    #                                 index_set=bus_attrs['names'],
    #                                 buses=buses,
    #                                 coordinate_type=CoordinateType.RECTANGULAR
    #                                 )

    ### declare angle difference limits on interconnected buses
    # libbranch.declare_ineq_angle_diff_branch_lbub(model=model,
    #                                               index_set=branch_attrs['names'],
    #                                               branches=branches,
    #                                               coordinate_type=CoordinateType.RECTANGULAR
    #                                               )

    libbranch.declare_socp_scw(model = model, index_set = branch_attrs['names'],
                                 branches = branches, 
                                 branch_attrs = branch_attrs)

    ### declare the generator cost objective
    libgen.declare_expression_pgqg_operating_cost(model=model,
                                                  index_set=gen_attrs['names'],
                                                  p_costs=gen_attrs['p_cost'],
                                                  q_costs=gen_attrs.get('q_cost', None)
                                                  )

    obj_expr = sum(model.pg_operating_cost[gen_name] for gen_name in model.pg_operating_cost)
    if include_feasibility_slack:
        obj_expr += penalty_expr
    if hasattr(model, 'qg_operating_cost'):
        obj_expr += sum(model.qg_operating_cost[gen_name] for gen_name in model.qg_operating_cost)

    model.obj = pe.Objective(expr=obj_expr)

    return model, md


if __name__ == '__main__':
    import os
    from egret.parsers.matpower_parser import create_ModelData
    from egret.common.solver_interface import _solve_model

    path = os.path.dirname(__file__)
    #filename = 'pglib_opf_case3_lmbd.m'
    #filename = 'pglib_opf_case5_pjm.m'
    #filename = 'pglib_opf_case14_ieee.m'
    filename = 'pglib_opf_case118_ieee.m'
    matpower_file = os.path.join(path, '../../downloads/pglib-opf-master/', filename)
    model_data = create_ModelData(matpower_file)
    print(len(dict(model_data.elements(element_type='generator')))) #will give the number of generators 
    #model,md = create_socp_acopf_model(model_data)
    #m, results = _solve_model(model,'ipopt',solver_tee = True)

    #model.pprint()



