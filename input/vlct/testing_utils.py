# This file defines some useful testing utilities shared between a couple of
# regression testing scripts

import os
import os.path
import subprocess

try:
  basestring
except NameError:
  basestring = str

import numpy as np

def prep_cur_dir(executable):
    cwd = os.getcwd()
    if cwd[-10:] == "input/vlct":
        os.chdir("../../")
        # the following is just for a possible error message
        orig_exec_path = os.path.join(cwd,"../../",executable)
    else:
        # the following is just for a possible error message
        orig_exec_path = os.path.join(cwd,executable)

    if not os.path.isfile(executable):
        raise RuntimeError("Can't locate the executable: " + orig_exec_path)

class _BaseCalcL1Norm(object):
    """
    Base class used to define configurable functor objecst that encapsulate the
    calls to the l1_error_norm script
    """
    def __init__(self, script_path, default_fields = None, verbose = False):
        self.script_path = script_path

        self._check_field_container(default_fields, "default_fields")
        self.default_fields = default_fields
        self._verbose = verbose

    def __call__(self, dir_name, *args, **kwargs):
        raise NotImplementedError("This must be implemented by subclasses")

    def _check_field_container(self, val, argname):
        """                                                                     
        Checks that a argument passed that holds a list of fields, is None, or
        it is a sequence of strings.
        """
        if val is None:
            return
        _err_prefix \
            = "{} must be None or a container of strings.".format(argname)
        if isinstance(val, basestring):
            raise ValueError(_err_prefix + "It cannot be a string.")
        if not all((isinstance(field,basestring) for field in val)):
            raise ValueError(_err_prefix + "It contains non-string entries.")

    def _comma_separated_fields(self, fields):
        """
        Returns a string of comma separated fields from the fields argument. If
        fields is None, then the list is determined from the default_fields.
        If both are None, then None is returned.
        """
        str_field_l = None
        if fields is not None:
            _check_field_container(self, fields, 'fields')
            str_field_l = ','.join(fields)
        elif self.default_fields is not None:
            str_field_l = ','.join(self.default_fields)
        return str_field_l

    def _run_command(self, command_list):
        """
        Responsible for actually running the command. 

        The elements of command_list should be the command to execute on the
        shell followed by a list of arguments.
        """
        command = ' '.join(command_list)
        if self._verbose:
           print("calling {:s}".format(command))
        return float(subprocess.check_output(command,shell=True))
      

class CalcSimL1Norm(_BaseCalcL1Norm):
    """
    Configurable functor object that encapsulates the call to the l1_error_norm
    script in simulation mode (where L1 norm is computed by comparing 2 
    simulation outputs)

    Parameters
    ----------
    script_path : str
        Path to the l1_error_norm script
    default_fields : list of strings, optional
        Default list of fields to use in L1 Error Norm calculation by default.
        This can been overwritten while calling the functor
    """
    def __init__(self,script_path, default_fields = None, verbose = False):
        super(CalcSimL1Norm, self).__init__(script_path, default_fields,
                                            verbose)

    def __call__(self, dir_name, dir_name2, res = None, fields = None):
        """
        Computes the norm of the L1 error vector for one snapshot with respect 
        to another snapshot

        Parameters
        ----------
        dir_name : str
            The path to the directory holding one of the twos snapshots
        dir_name2 : str
            The path to the directory holding the other snapshot
        res : int, optional
            If this is specified then the normalization of the norm of the l1 
            error vector is computed assuming that the shape is (res,res,res).
        fields : list of strings, optional
            If this is specified then these are the names of the fields 
            included in L1 error vector. If this is not specified, the fields
            listed in the default_fields attribute are used. If those are not
            specified either, then all the fields appearing in both snapshots 
            are used.
        """

        command_list = ["python", self.script_path, "sim", dir_name, dir_name2]

        if res is not None:
            if ((int(res) != res) or res<=0):
                raise ValueError("res must be a positive integer")
            command_list += ["-n", str(int(res))]

        str_field_l = self._comma_separated_fields(fields)
        if str_field_l is not None:
            command_list += ['-f',str_field_l]

        return self._run_command(command_list)

class CalcTableL1Norm(_BaseCalcL1Norm):
    """
    Configurable functor object that encapsulates the call to the l1_error_norm
    script in tablex, tabley, or tablez mode (L1 norm is computed by comparing
    a 1D slice of a simulation output along the specified axis to tabulated 
    values)

    Parameters
    ----------
    script_path : str
        Path to the l1_error_norm script
    default_fields : list of strings, optional
        Default list of fields to use in L1 Error Norm calculation by default.
        This can been overwritten while calling the functor
    default_ref_table : str, optional
        Default reference table to use in the comparisons
    verbose : bool, optional
        Default is False. When True, extra information is printed to stdout
    """
    def __init__(self, script_path, default_fields = None, 
                 default_ref_table = None, verbose = False):
        super(CalcTableL1Norm,self).__init__(script_path, default_fields,
                                             verbose)
        self.default_ref_table = default_ref_table

    def _handle_bkg_velocity(self,bkg_velocity, permute, axis):
        out = []
        if len(bkg_velocity) == 3:
            if permute and isinstance(bkg_velocity,list):
                if axis.lower() == 'y':
                    bkg_velocity = bkg_velocity[2:] + bkg_velocity[:2]
                elif axis == 'z':
                    bkg_velocity = bkg_velocity[1:] + bkg_velocity[:1]
            out.append('--bkg_velocity')
            out += [repr(float(elem)) for elem in bkg_velocity]
        elif len(bkg_velocity) != 0:
            raise ValueError("invoking --bkg_velocity, requires 3 values")
        return out

    def __call__(self, dir_name, axis, fields = None, ref_table = None,
                 std_dev = False, permute = False, reverse = False,
                 bkg_velocity = (), offset_soln = None):
        """
        Computes the L1 Error norms

        If bkg_velocity is a tuple, then the values are passed to the command 
        line as-is. The same holds if bkg_velocity is a list and permute is 
        False. However, if bkg_velocity is a list and permute is True, then 
        entries of bkg_velocity are permuted according to axis.
        """
        tab_modes = {'x':'tablex', 'y':'tabley', 'z':'tablez'}
        mode = tab_modes.get(axis.lower(),None)
        if mode is None:
            raise ValueError("The axis arg must be 'x','y',or 'z'")

        if ref_table is None:
            if self.default_ref_table is None:
                raise ValueError("No reference table was specified")
            ref_table = self.default_ref_table

        command_list = ["python", self.script_path, mode, ref_table, dir_name]

        str_field_l = self._comma_separated_fields(fields)
        if str_field_l is not None:
            command_list += ['-f',str_field_l]

        if permute:
            command_list.append('--permute')
        if std_dev:
            command_list.append('--std')
        if reverse:
            command_list.append('--reverse')

        command_list += self._handle_bkg_velocity(bkg_velocity, permute, axis)

        if offset_soln is not None:
            command_list.append('--offset_soln')
            command_list.append(repr(float(offset_soln)))

        return self._run_command(command_list)

def isclose(a,b,abs_tol = False):
    # just wraps numpy's isclose function
    # this is how much we allow things to be wrong (a little arbitrary)
    err = 1.e-13
    # it makes sense to allow for absolute errors in l1 norms since the
    # floating point differences will enter into the values that are being
    # subtracted
    if abs_tol:
        # slow wave requires slightly bigger tolerance (since it includes more
        # timesteps which allows the floating point errors to grow more)
        atol = 2.e-14
    else:
        atol = 0
    return np.isclose(a,b,rtol=err,atol=atol)

def standard_analyze(ref_val, L1Functor, target_path, test_case_name,
                     functor_kwargs = {}, target_path_args = (),
                     exact = False, verbose = False):
    """
    Computes the L1 Error norm and computes the result.

    This provides standardized output messages for L1 Error norms computed both
    by comparing values between snapshots in a simulation AND by comparing 
    snapshots to a tabulated solution.

    Parameters
    ----------
    ref_val : float or dict
        This is the floating point reference value.
    L1Functor: instance of CalcSimL1Norm or CalcTableL1Norm
        This functor wraps the call to the l1_error_norm script.
    test_case_name : string
        The name of the test case to report when printing the results
    target_path : string
        The path to the directory containing data to compute L1 Error Norm for.
    functor_kwargs : dict of keyword arguments, optional
        Optional kwargs to pass to L1Functor. Default is an empty dict
    target_path_args : tuple of arguments
        A list of arguments to pass to target_path.format. If this is empty, 
        (by default) then the format is method not called.
    exact : bool, optional
        Indicates whether the computed L1 Error norm should exactly match 
        ref_val. Default is False.
    verbose : bool, optional
        If False (default), just print messages denoting test failure
    """
    if len(target_path_args) != 0:
        target_path = target_path.format(*target_path_args)

    norm = L1Functor(dir_name = target_path, **functor_kwargs)

    std_dev = functor_kwargs.get('std_dev',False)
    if exact != std_dev:
        raise NotImplementedError()

    msg = None
    if exact and std_dev:
        passed = (norm == ref_val)
        if not passed:
            msg = ("FAILED: Standard Deviation of norm of L1 error of {case} "
                   "isn't correct\ncomputed: {norm} expected: {ref}")
        elif verbose:
            msg = ("PASSED: Standard Deviation of norm of L1 error of {case} "
                   "has the correct value of {norm}")
    elif (not exact):
        passed = isclose(norm, ref_val, abs_tol = True)
        if not passed:
            msg = ("FAILED: L1 error of {case} isn't correct\n"
                   "computed: {norm} expected: {ref}")
        elif verbose:
            msg = ("PASSED: L1 error of {case} within tolerance\n"
                   "computed: {norm} expected: {ref}")
    if msg is not None:
        print(msg.format(case=test_case_name, norm=repr(norm),
                         ref=repr(ref_val)))
    return True

class EnzoEWrapper:
    def __init__(self, executable, input_file_template):
        self.command_template = executable + ' ' + input_file_template
    def __call__(self, *args):
        if len(args) == 0:
            command = self.command_template
        else:
            command = self.command_template.format(*args)
        subprocess.call(command,shell=True)
