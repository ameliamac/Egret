#  ___________________________________________________________________________
#
#  EGRET: Electrical Grid Research and Engineering Tools
#  Copyright 2019 National Technology & Engineering Solutions of Sandia, LLC
#  (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S.
#  Government retains certain rights in this software.
#  This software is distributed under the Revised BSD License.
#  ___________________________________________________________________________

"""
This module contains the declarations for the modeling components
typically used for transmission lines
"""
import math
import pyomo.environ as pe
import egret.model_library.transmission.tx_calc as tx_calc
import egret.model_library.decl as decl
from egret.model_library.defn import FlowType, CoordinateType, ApproximationType, RelaxationType
from egret.data.model_data import zip_items


def declare_var_dva(model, index_set, **kwargs):
    """
    Create variable or the angle difference between interconnected bus pairs
    """
    decl.declare_var('dva', model=model, index_set=index_set, **kwargs)


def declare_var_pfl(model, index_set, **kwargs):
    """
    Create variable for the real part of the power loss in the transmission
    line
    """
    decl.declare_var('pfl', model=model, index_set=index_set, **kwargs)


def declare_var_pf(model, index_set, **kwargs):
    """
    Create variable for the real part of the power flow in the "from"
    end of the transmission line
    """
    decl.declare_var('pf', model=model, index_set=index_set, **kwargs)

def declare_expr_pf(model, index_set, **kwargs):
    """
    Create expression for the real part of the power flow in the "from"
    end of the transmission line
    """
    decl.declare_expr('pf', model=model, index_set=index_set, **kwargs)

def declare_var_qf(model, index_set, **kwargs):
    """
    Create variable for the imaginary part of the power flow in the "from"
    end of the transmission line
    """
    decl.declare_var('qf', model=model, index_set=index_set, **kwargs)


def declare_var_pt(model, index_set, **kwargs):
    """
    Create variable for the real part of the power flow in the "to"
    end of the transmission line
    """
    decl.declare_var('pt', model=model, index_set=index_set, **kwargs)


def declare_var_qt(model, index_set, **kwargs):
    """
    Create variable for the imaginary part of the power flow in the "to"
    end of the transmission line
    """
    decl.declare_var('qt', model=model, index_set=index_set, **kwargs)


def declare_var_ifr(model, index_set, **kwargs):
    """
    Create variable for the real part of the current flow in the "from"
    end of the transmission line
    """
    decl.declare_var('ifr', model=model, index_set=index_set, **kwargs)


def declare_var_ifj(model, index_set, **kwargs):
    """
    Create variable for the imaginary part of the current flow in the "from"
    end of the transmission line
    """
    decl.declare_var('ifj', model=model, index_set=index_set, **kwargs)


def declare_var_itr(model, index_set, **kwargs):
    """
    Create variable for the real part of the current flow in the "to"
    end of the transmission line
    """
    decl.declare_var('itr', model=model, index_set=index_set, **kwargs)


def declare_var_itj(model, index_set, **kwargs):
    """
    Create variable for the imaginary part of the current flow in the "to"
    end of the transmission line
    """
    decl.declare_var('itj', model=model, index_set=index_set, **kwargs)


def declare_eq_branch_dva(model, index_set, branches):
    """
    Create the equality constraints for the angle difference
    in the branch
    """
    m = model

    con_set = decl.declare_set("_con_eq_branch_dva_set", model, index_set)

    m.eq_dva_branch = pe.Constraint(con_set)
    for branch_name in con_set:
        branch = branches[branch_name]

        from_bus = branch['from_bus']
        to_bus = branch['to_bus']

        shift = 0.0
        if branch['branch_type'] == 'transformer':
            shift = math.radians(branch['transformer_phase_shift'])

        m.eq_dva_branch[branch_name] = \
            m.dva[branch_name] == \
            m.va[from_bus] - m.va[to_bus] + shift


def declare_expr_c(model, index_set, coordinate_type=CoordinateType.POLAR):
    """
    Create expression for the nonlinear, nonconvex term based on cosine
    of the phase angle difference (polar) or bilinear voltages (rectangular)
    """
    m = model
    expr_set = decl.declare_set('_expr_c', model, index_set)
    m.c = pe.Expression(expr_set)

    if coordinate_type == CoordinateType.RECTANGULAR:
        for from_bus, to_bus in expr_set:
            m.c[(from_bus,to_bus)] = m.vr[from_bus]*m.vr[to_bus] + m.vj[from_bus]*m.vj[to_bus]
    elif coordinate_type == CoordinateType.POLAR:
        for from_bus, to_bus in expr_set:
            m.c[(from_bus,to_bus)] = m.vm[from_bus]*m.vm[to_bus]*pe.cos(m.va[from_bus]-m.va[to_bus])


def declare_expr_s(model, index_set, coordinate_type=CoordinateType.POLAR):
    """
    Create expression for the nonlinear, nonconvex term based on cosine
    of the phase angle difference (polar) or bilinear voltages (rectangular)
    """
    m = model
    expr_set = decl.declare_set('_expr_s', model, index_set)
    m.s = pe.Expression(expr_set)


    if coordinate_type == CoordinateType.RECTANGULAR:
        for from_bus, to_bus in expr_set:
            m.s[(from_bus,to_bus)] = m.vj[from_bus]*m.vr[to_bus] - m.vr[from_bus]*m.vj[to_bus]
    elif coordinate_type == CoordinateType.POLAR:
        for from_bus, to_bus in expr_set:
            m.s[(from_bus,to_bus)] = m.vm[from_bus]*m.vm[to_bus]*pe.sin(m.va[from_bus]-m.va[to_bus])


def declare_var_c(model, index_set, **kwargs):
    """
    Create variable for the for the bilinear voltages in SOCP formulation
    """
    # if user provides bounds as dict of tuple, translate it
    # into something that Pyomo understands
    varname = 'c'
    if kwargs and 'bounds' in kwargs and isinstance(kwargs['bounds'], dict):
        d = kwargs['bounds']
        def _c_bounds(model,j,k):
            return(d[j,k]) #var c has tuple keys
        kwargs['bounds'] = _c_bounds

    # create var if index set is None
    if index_set is None:
        model.add_component(varname, pe.Var(**kwargs))
    # transform the index set into a Pyomo Set
    else:
        pyomo_index_set = pe.Set(initialize=index_set, ordered=True)
        model.add_component("_var_{}_index_set".format(varname), pyomo_index_set)

        # now create the var
        model.add_component(varname, pe.Var(pyomo_index_set, **kwargs))

def declare_var_s(model, index_set, **kwargs):
    """
    Create variable for the for the bilinear voltages in SOCP formulation
    """
    varname = 's'
    if kwargs and 'bounds' in kwargs and isinstance(kwargs['bounds'], dict):
        d = kwargs['bounds']
        def _c_bounds(model,j,k):
            return(d[j,k]) #var s has tuple keys
        kwargs['bounds'] = _c_bounds

    # create var if index set is None
    if index_set is None:
        model.add_component(varname, pe.Var(**kwargs))
    # transform the index set into a Pyomo Set
    else:
        pyomo_index_set = pe.Set(initialize=index_set, ordered=True)
        model.add_component("_var_{}_index_set".format(varname), pyomo_index_set)

        # now create the var
        model.add_component(varname, pe.Var(pyomo_index_set, **kwargs))



def declare_eq_branch_current(model, index_set, branches, coordinate_type=CoordinateType.RECTANGULAR):
    """
    Create the equality constraints for the real and imaginary current
    in the branch
    """
    assert(coordinate_type != CoordinateType.POLAR
           and "Branch current in polar coordinates not implemented.")

    m = model
    con_set = decl.declare_set("_con_eq_branch_current_set", model, index_set)

    m.eq_ifr_branch = pe.Constraint(con_set)
    m.eq_ifj_branch = pe.Constraint(con_set)
    m.eq_itr_branch = pe.Constraint(con_set)
    m.eq_itj_branch = pe.Constraint(con_set)
    for branch_name in con_set:
        branch = branches[branch_name]

        from_bus = branch['from_bus']
        to_bus = branch['to_bus']

        g = tx_calc.calculate_conductance(branch)
        b = tx_calc.calculate_susceptance(branch)
        bc = branch['charging_susceptance']
        tau = 1.0
        shift = 0.0

        if branch['branch_type'] == 'transformer':
            tau = branch['transformer_tap_ratio']
            shift = math.radians(branch['transformer_phase_shift'])

        g11 = g / tau**2
        g12 = (g * math.cos(shift) - b * math.sin(shift)) / tau
        g21 = (g * math.cos(shift) + b * math.sin(shift)) / tau
        g22 = g

        b11 = (b + bc / 2) / tau**2
        b12 = (b * math.cos(shift) + g*math.sin(shift)) / tau
        b21 = (b * math.cos(shift) - g*math.sin(shift)) / tau
        b22 = b + bc / 2

        m.eq_ifr_branch[branch_name] = \
            m.ifr[branch_name] == \
            g11 * m.vr[from_bus] - g12 * m.vr[to_bus] - (b11 * m.vj[from_bus] - b12 * m.vj[to_bus])

        m.eq_ifj_branch[branch_name] = \
            m.ifj[branch_name] == \
            g11 * m.vj[from_bus] - g12 * m.vj[to_bus] + (b11 * m.vr[from_bus] - b12 * m.vr[to_bus])

        m.eq_itr_branch[branch_name] = \
            m.itr[branch_name] == \
            -(g21 * m.vr[from_bus] - g22 * m.vr[to_bus] - (b21 * m.vj[from_bus] - b22 * m.vj[to_bus]))

        m.eq_itj_branch[branch_name] = \
            m.itj[branch_name] == \
            -(g21 * m.vj[from_bus] - g22 * m.vj[to_bus] + (b21 * m.vr[from_bus] - b22 * m.vr[to_bus]))


def declare_eq_branch_power(model, index_set, branches, branch_attrs, coordinate_type=CoordinateType.POLAR):
    """
    Create the equality constraints for the real and reactive power
    in the branch
    """
    m = model

    bus_pairs = zip_items(branch_attrs['from_bus'],branch_attrs['to_bus'])
    unique_bus_pairs = list(set([val for idx,val in bus_pairs.items()]))
    declare_expr_c(model,unique_bus_pairs,coordinate_type)
    declare_expr_s(model,unique_bus_pairs,coordinate_type)

    con_set = decl.declare_set("_con_eq_branch_power_set", model, index_set)

    m.eq_pf_branch = pe.Constraint(con_set)
    m.eq_pt_branch = pe.Constraint(con_set)
    m.eq_qf_branch = pe.Constraint(con_set)
    m.eq_qt_branch = pe.Constraint(con_set)
    for branch_name in con_set:
        branch = branches[branch_name]

        from_bus = branch['from_bus']
        to_bus = branch['to_bus']

        if coordinate_type == CoordinateType.POLAR:
            vmsq_from_bus = m.vm[from_bus]**2
            vmsq_to_bus = m.vm[to_bus] ** 2
        elif coordinate_type == CoordinateType.RECTANGULAR:
            vmsq_from_bus = m.vr[from_bus]**2 + m.vj[from_bus]**2
            vmsq_to_bus = m.vr[to_bus] ** 2 + m.vj[to_bus] ** 2

        g = tx_calc.calculate_conductance(branch)
        b = tx_calc.calculate_susceptance(branch)
        bc = branch['charging_susceptance']
        tau = 1.0
        shift = 0.0

        if branch['branch_type'] == 'transformer':
            tau = branch['transformer_tap_ratio']
            shift = math.radians(branch['transformer_phase_shift'])

        g11 = g / tau ** 2
        g12 = g * math.cos(shift) / tau
        g21 = g * math.sin(shift) / tau
        g22 = g

        b11 = (b + bc / 2) / tau ** 2
        b12 = b * math.cos(shift) / tau
        b21 = b * math.sin(shift) / tau
        b22 = b + bc / 2

        m.eq_pf_branch[branch_name] = \
            m.pf[branch_name] == \
            g11 * vmsq_from_bus - \
            (g12 * m.c[(from_bus,to_bus)] +
             g21 * m.s[(from_bus,to_bus)] +
             b12 * m.s[(from_bus,to_bus)] -
             b21 * m.c[(from_bus,to_bus)])

        m.eq_pt_branch[branch_name] = \
            m.pt[branch_name] == \
            g22 * vmsq_to_bus - \
            (g12 * m.c[(from_bus,to_bus)] +
             g21 * m.s[(from_bus,to_bus)] -
             b12 * m.s[(from_bus,to_bus)] +
             b21 * m.c[(from_bus,to_bus)])

        m.eq_qf_branch[branch_name] = \
            m.qf[branch_name] == \
            -b11 * vmsq_from_bus + \
            (b12 * m.c[(from_bus,to_bus)] +
             b21 * m.s[(from_bus,to_bus)] -
             g12 * m.s[(from_bus,to_bus)] +
             g21 * m.c[(from_bus,to_bus)])

        m.eq_qt_branch[branch_name] = \
            m.qt[branch_name] == \
            -b22 * vmsq_to_bus + \
            (b12 * m.c[(from_bus,to_bus)] +
             b21 * m.s[(from_bus,to_bus)] +
             g12 * m.s[(from_bus,to_bus)] -
             g21 * m.c[(from_bus,to_bus)])


def declare_eq_branch_power_btheta_approx(model, index_set, branches, approximation_type=ApproximationType.BTHETA):
    """
    Create the equality constraints for power (from BTHETA approximation)
    in the branch
    """
    m = model

    con_set = decl.declare_set("_con_eq_branch_power_btheta_approx_set", model, index_set)

    m.eq_pf_branch = pe.Constraint(con_set)
    for branch_name in con_set:
        branch = branches[branch_name]

        from_bus = branch['from_bus']
        to_bus = branch['to_bus']

        tau = 1.0
        shift = 0.0
        if branch['branch_type'] == 'transformer':
            tau = branch['transformer_tap_ratio']
            shift = math.radians(branch['transformer_phase_shift'])

        if approximation_type == ApproximationType.BTHETA:
            x = branch['reactance']
            b = -1/(tau*x)
        elif approximation_type == ApproximationType.BTHETA_LOSSES:
            b = tx_calc.calculate_susceptance(branch)/tau

        m.eq_pf_branch[branch_name] = \
            m.pf[branch_name] == \
            b * (m.va[from_bus] - m.va[to_bus] + shift)


def declare_eq_branch_loss_btheta_approx(model, index_set, branches, relaxation_type = RelaxationType.NONE):
    """
    Create the equality constraints for losses (from BTHETA approximation)
    in the branch
    """
    m = model

    con_set = decl.declare_set("_con_eq_branch_loss_btheta_approx_set", model, index_set)

    m.eq_pfl_branch = pe.Constraint(con_set)
    for branch_name in con_set:
        branch = branches[branch_name]

        tau = 1.0
        if branch['branch_type'] == 'transformer':
            tau = branch['transformer_tap_ratio']
        g = tx_calc.calculate_conductance(branch)/tau

        if relaxation_type == RelaxationType.NONE:
            m.eq_pfl_branch[branch_name] = \
                m.pfl[branch_name] == \
                g * (m.dva[branch_name])**2
        elif relaxation_type == RelaxationType.SOC:
            m.eq_pfl_branch[branch_name] = \
                m.pfl[branch_name] >= \
                g * (m.dva[branch_name])**2


def get_power_flow_expr_ptdf_approx(model, branch_name, PTDF, rel_ptdf_tol=None, abs_ptdf_tol=None):
    """
    Create a pyomo power flow expression from PTDF matrix
    """

    if rel_ptdf_tol is None:
        rel_ptdf_tol = 0.
    if abs_ptdf_tol is None:
        abs_ptdf_tol = 0.

    expr = PTDF.get_branch_phase_shift(branch_name)

    max_coef = PTDF.get_branch_ptdf_abs_max(branch_name)
    ## This case is weird, but could happen, causing divison by 0 below
    if max_coef == 0:
        return expr
    ## NOTE: It would be easy to hold on to the 'ptdf' dictionary here,
    ##       if we wanted to
    for bus_name, coef in PTDF.get_branch_ptdf_iterator(branch_name):
        if abs(coef) < abs_ptdf_tol:
            ## no point in excuting the rest of the for loop
            continue
        if abs(coef)/max_coef < rel_ptdf_tol:
            ## no point in excuting the rest of the for loop
            continue
        phi_adjust = PTDF.get_bus_phi_adj(bus_name)

        expr += coef*model.p_nw[bus_name]+coef*phi_adjust

    return expr


def declare_eq_branch_power_ptdf_approx(model, index_set, PTDF, rel_ptdf_tol=None, abs_ptdf_tol=None):
    """
    Create the equality constraints or expressions for power (from PTDF 
    approximation) in the branch
    """

    m = model

    con_set = decl.declare_set("_con_eq_branch_power_ptdf_approx_set", model, index_set)

    pf_is_var = isinstance(m.pf, pe.Var)

    if pf_is_var:
        m.eq_pf_branch = pe.Constraint(con_set)
    else:
        if not isinstance(m.pf, pe.Expression):
            raise Exception("Unrecognized type for m.pf", m.pf.pprint())

    for branch_name in con_set:
        expr = \
            get_power_flow_expr_ptdf_approx(m, branch_name, PTDF, rel_ptdf_tol=rel_ptdf_tol, abs_ptdf_tol=abs_ptdf_tol)

        if pf_is_var:
            m.eq_pf_branch[branch_name] = \
                m.pf[branch_name] == expr
        else:
            m.pf[branch_name] = expr


def get_branch_loss_expr_ptdf_approx(model, branch_name, PTDF, rel_ptdf_tol=None, abs_ptdf_tol=None): 
    """
    Create a pyomo power flow loss expression from PTDF matrix
    """
    if rel_ptdf_tol is None:
        rel_ptdf_tol = 0.
    if abs_ptdf_tol is None:
        abs_ptdf_tol = 0.

    expr = PTDF.get_branch_losses_phase_shift(branch_name)
    expr += PTDF.get_branch_ldf_c(branch_name)

    max_coef = PTDF.get_branch_ldf_abs_max(branch_name)

    if max_coef == 0:
        return expr

    ## NOTE: It would be easy to hold on to the 'ldf' dictionary here,
    ##       if we wanted to
    for bus_name, coef in PTDF.get_branch_ldf_iterator(branch_name):
        if abs(coef) < abs_ptdf_tol:
            ## no point in excuting the rest of the for loop
            continue
        if abs(coef)/max_coef < rel_ptdf_tol:
            ## no point in excuting the rest of the for loop
            continue
        phi_losses_adjust = PTDF.get_bus_phi_losses_adj(bus_name)

        expr += coef*model.p_nw[bus_name]+coef*phi_losses_adjust

    return expr


def declare_eq_branch_loss_ptdf_approx(model, index_set, PTDF, rel_ptdf_tol=None, abs_ptdf_tol=None):
    """
    Create the equality constraints or expressions for losses (from PTDF 
    approximation) in the branch
    """
    m = model

    con_set = decl.declare_set("_con_eq_branch_loss_ptdf_approx_set", model, index_set)
    pfl_is_var = isinstance(m.pfl, pe.Var)
    if pfl_is_var:
        m.eq_pfl_branch = pe.Constraint(con_set)
    else:
        if not isinstance(m.pfl, pe.Expression):
            raise Exception("Unrecognized type for m.pfl", m.pfl.pprint())

    for branch_name in con_set:
        expr = \
            get_branch_loss_expr_ptdf_approx(m, branch_name, PTDF, rel_ptdf_tol=rel_ptdf_tol, abs_ptdf_tol=abs_ptdf_tol)

        if pfl_is_var:
            m.eq_pfl_branch[branch_name] = \
                m.pfl[branch_name] == expr
        else:
            m.pfl[branch_name] = expr

def declare_ineq_s_branch_thermal_limit(model, index_set,
                                        branches, s_thermal_limits,
                                        flow_type=FlowType.POWER):
    """
    Create the inequality constraints for the branch thermal limits
    based on the power variables.
    """
    m = model
    con_set = decl.declare_set('_con_ineq_s_branch_thermal_limit',
                               model=model, index_set=index_set)

    m.ineq_sf_branch_thermal_limit = pe.Constraint(con_set)
    m.ineq_st_branch_thermal_limit = pe.Constraint(con_set)

    if flow_type == FlowType.CURRENT:
        for branch_name in con_set:
            if s_thermal_limits[branch_name] is None:
                continue

            from_bus = branches[branch_name]['from_bus']
            to_bus = branches[branch_name]['to_bus']
            m.ineq_sf_branch_thermal_limit[branch_name] = \
                (m.vr[from_bus] ** 2 + m.vj[from_bus] ** 2) * (m.ifr[branch_name] ** 2 + m.ifj[branch_name] ** 2) \
                <= s_thermal_limits[branch_name] ** 2
            m.ineq_st_branch_thermal_limit[branch_name] = \
                (m.vr[to_bus] ** 2 + m.vj[to_bus] ** 2) * (m.itr[branch_name] ** 2 + m.itj[branch_name] ** 2) \
                <= s_thermal_limits[branch_name] ** 2
    elif flow_type == FlowType.POWER:
        for branch_name in con_set:
            if s_thermal_limits[branch_name] is None:
                continue

            m.ineq_sf_branch_thermal_limit[branch_name] = \
                m.pf[branch_name] ** 2 + m.qf[branch_name] ** 2 \
                <= s_thermal_limits[branch_name] ** 2
            m.ineq_st_branch_thermal_limit[branch_name] = \
                m.pt[branch_name] ** 2 + m.qt[branch_name] ** 2 \
                <= s_thermal_limits[branch_name] ** 2


def declare_ineq_p_branch_thermal_lbub(model, index_set,
                                        branches, p_thermal_limits,
                                        approximation_type=ApproximationType.BTHETA):
    """
    Create the inequality constraints for the branch thermal limits
    based on the power variables or expressions.
    """
    m = model
    con_set = decl.declare_set('_con_ineq_p_branch_thermal_lbub',
                               model=model, index_set=index_set)

    m.ineq_pf_branch_thermal_lb = pe.Constraint(con_set)
    m.ineq_pf_branch_thermal_ub = pe.Constraint(con_set)

    if approximation_type == ApproximationType.BTHETA or \
            (approximation_type == ApproximationType.PTDF and \
            isinstance(m.pf, pe.Expression)):
        for branch_name in con_set:
            if p_thermal_limits[branch_name] is None:
                continue

            m.ineq_pf_branch_thermal_lb[branch_name] = \
                -p_thermal_limits[branch_name] <= m.pf[branch_name]

            m.ineq_pf_branch_thermal_ub[branch_name] = \
                m.pf[branch_name] <= p_thermal_limits[branch_name]


def declare_ineq_angle_diff_branch_lbub(model, index_set,
                                        branches,
                                        coordinate_type=CoordinateType.POLAR):
    """
    Create the inequality constraints for the angle difference
    bounds between interconnected buses.
    """
    m = model
    con_set = decl.declare_set('_con_ineq_angle_diff_branch_lbub',
                               model=model, index_set=index_set)

    m.ineq_angle_diff_branch_lb = pe.Constraint(con_set)
    m.ineq_angle_diff_branch_ub = pe.Constraint(con_set)

    if coordinate_type == CoordinateType.POLAR:
        for branch_name in con_set:
            from_bus = branches[branch_name]['from_bus']
            to_bus = branches[branch_name]['to_bus']

            m.ineq_angle_diff_branch_lb[branch_name] = \
                branches[branch_name]['angle_diff_min'] <= m.va[from_bus] - m.va[to_bus]
            m.ineq_angle_diff_branch_ub[branch_name] = \
                m.va[from_bus] - m.va[to_bus] <= branches[branch_name]['angle_diff_max']
    elif coordinate_type == CoordinateType.RECTANGULAR:
        for branch_name in con_set:
            from_bus = branches[branch_name]['from_bus']
            to_bus = branches[branch_name]['to_bus']

            m.ineq_angle_diff_branch_lb[branch_name] = \
                branches[branch_name]['angle_diff_min'] <= pe.atan(m.vj[from_bus]/m.vr[from_bus]) \
                - pe.atan(m.vj[to_bus]/m.vr[to_bus])
            m.ineq_angle_diff_branch_ub[branch_name] = \
                pe.atan(m.vj[from_bus] / m.vr[from_bus]) \
                - pe.atan(m.vj[to_bus] / m.vr[to_bus]) <= branches[branch_name]['angle_diff_max']


def declare_eq_branch_power_socp(model, index_set, branches, branch_attrs, coordinate_type=CoordinateType.POLAR):
    """
    Create the equality constraints for the real and reactive power
    in the branch with SOCP variables s,c,w
    """
    m = model

    bus_pairs = zip_items(branch_attrs['from_bus'],branch_attrs['to_bus'])
    unique_bus_pairs = list(set([val for idx,val in bus_pairs.items()]))
    #exprs are not needed because we have variables for these now
    # declare_expr_c(model,unique_bus_pairs,coordinate_type)
    # declare_expr_s(model,unique_bus_pairs,coordinate_type)

    con_set = decl.declare_set("_con_eq_branch_power_set", model, index_set)

    m.eq_pf_branch = pe.Constraint(con_set)
    m.eq_pt_branch = pe.Constraint(con_set)
    m.eq_qf_branch = pe.Constraint(con_set)
    m.eq_qt_branch = pe.Constraint(con_set)
    for branch_name in con_set:
        branch = branches[branch_name]

        from_bus = branch['from_bus']
        to_bus = branch['to_bus']

        #These are now represented by w
        # if coordinate_type == CoordinateType.POLAR:
        #     vmsq_from_bus = m.vm[from_bus]**2
        #     vmsq_to_bus = m.vm[to_bus] ** 2
        # elif coordinate_type == CoordinateType.RECTANGULAR:
        #     vmsq_from_bus = m.vr[from_bus]**2 + m.vj[from_bus]**2
        #     vmsq_to_bus = m.vr[to_bus] ** 2 + m.vj[to_bus] ** 2

        g = tx_calc.calculate_conductance(branch)
        b = tx_calc.calculate_susceptance(branch)
        bc = branch['charging_susceptance']
        tau = 1.0
        shift = 0.0

        if branch['branch_type'] == 'transformer':
            tau = branch['transformer_tap_ratio']
            shift = math.radians(branch['transformer_phase_shift'])

        g11 = g / tau ** 2
        g12 = g * math.cos(shift) / tau
        g21 = g * math.sin(shift) / tau
        g22 = g

        b11 = (b + bc / 2) / tau ** 2
        b12 = b * math.cos(shift) / tau
        b21 = b * math.sin(shift) / tau
        b22 = b + bc / 2

        m.eq_pf_branch[branch_name] = \
            m.pf[branch_name] == \
            g11 * m.w[from_bus] - \
            (g12 * m.c[(from_bus,to_bus)] +
             g21 * m.s[(from_bus,to_bus)] +
             b12 * m.s[(from_bus,to_bus)] -
             b21 * m.c[(from_bus,to_bus)])

        m.eq_pt_branch[branch_name] = \
            m.pt[branch_name] == \
            g22 * m.w[to_bus] - \
            (g12 * m.c[(from_bus,to_bus)] +
             g21 * m.s[(from_bus,to_bus)] -
             b12 * m.s[(from_bus,to_bus)] +
             b21 * m.c[(from_bus,to_bus)])

        m.eq_qf_branch[branch_name] = \
            m.qf[branch_name] == \
            -b11 * m.w[from_bus] + \
            (b12 * m.c[(from_bus,to_bus)] +
             b21 * m.s[(from_bus,to_bus)] -
             g12 * m.s[(from_bus,to_bus)] +
             g21 * m.c[(from_bus,to_bus)])

        m.eq_qt_branch[branch_name] = \
            m.qt[branch_name] == \
            -b22 * m.w[to_bus] + \
            (b12 * m.c[(from_bus,to_bus)] +
             b21 * m.s[(from_bus,to_bus)] +
             g12 * m.s[(from_bus,to_bus)] -
             g21 * m.c[(from_bus,to_bus)])


def declare_socp_scw(model, index_set, branches, branch_attrs):
    """
    Create the SOCP constraints c_bk^2+s_bk^2<= w_bb*w_kk
    """
    m = model

    con_set = decl.declare_set("_con_socp_scw_set", model, index_set)
    m.socp_scw = pe.Constraint(con_set)

    bus_pairs = zip_items(branch_attrs['from_bus'],branch_attrs['to_bus'])
    unique_bus_pairs = list(set([val for idx,val in bus_pairs.items()]))

    for branch_name in con_set:
        branch = branches[branch_name]
        from_bus = branch['from_bus']
        to_bus = branch['to_bus']

        m.socp_scw[branch_name] = \
        m.c[(from_bus,to_bus)]**2 + m.s[(from_bus,to_bus)]**2 <= \
        m.w[from_bus] * m.w[to_bus]





