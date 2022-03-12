import argparse
import io
import json

import numpy as np
import yt

from l1_error_norm import StoreCommaListAction
from test_report import create_test_report, create_dummy_report

_TABLE_DTYPE_ARGS = [
    ('name',  object),
    ('min', 'f8'),
    ('min_xloc', 'f8'),
    ('min_yloc', 'f8'),
    ('min_zloc', 'f8'),
    ('max', 'f8'),
    ('max_xloc', 'f8'),
    ('max_yloc', 'f8'),
    ('max_zloc', 'f8'),
    ('mean', 'f8'),
    ('variance', 'f8')
]

_FIXED_COL_UNITS = dict((col, 'code_length') for col,_ in _TABLE_DTYPE_ARGS \
                        if col[-5:] in ['_xloc', '_yloc', '_zloc'])

# sequence of tuples, mapping a yt derived_quantities to table columns.
_QUANTITY_ENTRY_SETS = [
    ('min_location', ('min', 'min_xloc', 'min_yloc', 'min_zloc'), {}),
    ('max_location', ('max', 'max_xloc', 'max_yloc', 'max_zloc'), {}),
    ('weighted_variance', ('variance', 'mean'),
     {'weight' : ("gas", "cell_volume")}),
]

_FMT_DICT = {object : 's', 'f8' : '.15e'}
_TABLE_FMT = ','.join('%'+_FMT_DICT[dtype] for _, dtype in _TABLE_DTYPE_ARGS)
_TABLE_COLUMNS = [col_name for col_name,_ in _TABLE_DTYPE_ARGS]
_TABLE_DTYPE = np.dtype(_TABLE_DTYPE_ARGS)

def construct_new_field_table(field_names):
    array = np.zeros((len(field_names),), dtype = _TABLE_DTYPE)
    for i,name in enumerate(field_names):
        array[i]['name'] = name
    return array

def _get_sim_props(ds, field_names):
    # determine the code unit definitions
    code_unit_defs = {}
    for unit in ds.unit_registry.keys():
        if 'code_' != unit[:5]:
            continue
        converted = ds.quan(1.0,unit).in_base('cgs')
        code_unit_defs[unit] = [float(converted.v), str(converted.units)]

    # determine the default units of each field
    fluid_type = None
    for elem in ds.fluid_types:
        if elem in ['enzoe','enzop']:
            assert fluid_type is None
            fluid_type = elem
    assert fluid_type is not None

    ds.field_list # this is necessary for initializing ds.field_info
    field_units = {}
    for field in field_names:
        dflt_units = ds.field_info[(fluid_type,field)].units
        # it may not be necessary to complete the following line
        converted = str(ds.quan(1.0, dflt_units).in_base('code').units)
        field_units[field] = converted
    return {'code unit definitions' : code_unit_defs,
            'field units' : field_units}

def measure_field_summary(target_path, field_names):
    target_ds = yt.load(target_path)
    sim_props = _get_sim_props(target_ds, field_names)

    computed_cols = set() # initialized for a sanity check!

    # initialize columns
    field_table = construct_new_field_table(field_names)
    computed_cols.add('name') # for sanity check

    # now, let's compute the field properties
    ad = target_ds.all_data()
    for derived_quantity, colnames, kwargs in _QUANTITY_ENTRY_SETS:
        for colname in colnames:
            assert colname not in computed_cols # sanity check
            computed_cols.add(colname)

        func = getattr(ad.quantities,derived_quantity)
        for row_ind, field_name in enumerate(field_names):
            rslt = func(field_name, **kwargs)
            for i,val in enumerate(rslt):
                colname = colnames[i]
                if colnames in _FIXED_COL_UNITS:
                    normalized = val.to(_FIXED_COL_UNITS[colnames])
                else:
                    normalized = val.in_base('code').v
                field_table[row_ind][colname] = normalized

    expected_cols = field_table.dtype.fields.keys()
    assert len(computed_cols.symmetric_difference(expected_cols)) == 0

    return field_table,sim_props

def write_field_summary(fname, field_table, sim_props = None):
    if hasattr(fname, 'write'):
        f = fname
        close_file = False
    else:
        f = open(fname, 'w')
        close_file = True

    if sim_props is not None:
        f.write('# ')
        f.write(json.dumps(sim_props))
        f.write('\n#\n')

    np.savetxt(f, field_table, fmt = _TABLE_FMT,
               header = ','.join(_TABLE_COLUMNS))

    if close_file:
        f.close()

def read_field_summary(fname):
    str_converter = lambda arg : arg.decode("utf-8")
    return np.genfromtxt(fname, dtype = _TABLE_DTYPE,
                         converters = {0: str_converter}, delimiter = ',')

def _check_consistent_cols_and_names(cur_field_table, ref_field_table, tr):
    # check for consistent columns
    cur_cols = cur_field_table.dtype.fields.keys()
    ref_cols = ref_field_table.dtype.fields.keys()
    consistent_cols = len(set(cur_cols).symmetric_difference(ref_cols)) == 0
    if not consistent_cols:
        tr.fail(
            ('Column Mismatch. Current table has the columns {!r} and the '
             'reference table has {!r}').format(list(cur_cols), list(ref_cols))
        )


    # check for consistent field names (for now, report different field
    # orderings as a problem)
    cur_names = cur_field_table['name']
    ref_names = ref_field_table['name']
    consistent_fields = ((len(cur_field_table) == len(ref_field_table)) and
                         (cur_names == ref_names).all())
    if not consistent_fields:
        tr.fail(
            ('Field Name Mismatch; current table stores data for {!r} while '
             'the reference stores data for {!r}').format(cur_names.tolist(),
                                                          ref_names.tolist())
        )
    return consistent_cols and consistent_fields

def test_equivalent_field_tables(cur_field_table, ref_field_table,
                                 atol = 1e-15, rtol = 1e-14,
                                 test_report = None):
    if test_report is None:
        tr = DummyReport()
    else:
        tr = test_report

    tr.write(
        'Comparing Field Tables (atol = {}, rtol = {})\n'.format(atol,rtol)
    )

    # first check that table properties are consistent
    consistent_table_prop = _check_consistent_cols_and_names(cur_field_table,
                                                             ref_field_table,
                                                             tr)
    if not consistent_table_prop: # exit early when inconsistent
        tr.incomplete('Aborting further comparisons')
        return False

    tr.write(
        'General table properties are consistent. Proceeding with comparison\n'
    )

    # Now actually verify equality of properties
    colnames = tuple(ref_field_table.dtype.fields.keys())
    comparison_arr = np.empty((len(cur_field_table), len(colnames)),
                              dtype = np.bool)
    for j,col in enumerate(colnames):
        if col == "name":
            comparison_arr[:,j] = True
        else:
            cur_vals = cur_field_table[col]
            ref_vals = ref_field_table[col]
            comparison_arr[:,j] = np.isclose(cur_vals, ref_vals, rtol = rtol,
                                             atol = atol, equal_nan = False)

    all_consistent_vals = True
    # Finally check comparison results on a row-by-row basis
    for i, name in enumerate(cur_field_table['name']):
        comparison_passes = comparison_arr[i,:].all()
        if comparison_passes:
            tr.passing(
                'All Summary properties for "{}" are consistent'.format(name)
            )
        else:
            num_consistent_props = 0
            inconsistent_comparisons = []
            for j,col in enumerate(colnames):
                if col == "name":
                    continue
                elif comparison_arr[i,j]:
                    num_consistent_props += 1
                    continue
                inconsistent_comparisons.append(
                    "'{}':\t cur = {:.15e}\tref = {:.15e}".format(
                        col, cur_field_table[i][col], ref_field_table[i][col]
                    )
                )
            num_props = num_consistent_props + len(inconsistent_comparisons)
            tmp='{} of {} summary properties are inconsistent for "{}"'.format(
                num_props - num_consistent_props, num_props, name
            )
            tr.fail('\n    '.join([tmp] + inconsistent_comparisons))
            all_consistent_vals = False
    return all_consistent_vals

def _main_measure(args):
    # Program to use by the measure subcommand. Just construct and report the
    # field summary table
    print("Measuring the Field Summary Properties")
    field_table, sim_props = measure_field_summary(
        args.target_path, field_names = args.fields
    )

    if args.output == '-':
        output = io.StringIO()
        sim_props = None
    else:
        print("Writing Field Summary Table to " + args.output)
        output = open(args.output, 'w')

    write_field_summary(output, field_table, sim_props = sim_props)

    if isinstance(output, io.StringIO):
        print("Printing Field Summary Table:\n")
        print(output.getvalue())
    output.close()


def _main_cmp(args):
    # Program to use by the cmp subcommand. Just construct the field summary
    # table and compare against a reference
    print("Loading the reference table to identify the fields that are to be "
          "summarized")
    ref_field_table = read_field_summary(args.ref)

    field_names = ref_field_table['name'].tolist()
    print("Measuring the Field Summary Properties")
    cur_field_table, sim_props = measure_field_summary(
        args.target_path, field_names = field_names
    )

    if args.report is None:
        report_creator = create_dummy_report
    else:
        report_creator = create_test_report

    with report_creator(args.report, clobber = True) as tr:
        print("Comparing field summary tables")
        test_rslt = test_equivalent_field_tables(cur_field_table,
                                                 ref_field_table,
                                                 atol = args.atol,
                                                 rtol = args.rtol,
                                                 test_report = tr)
        if test_rslt:
            print("Field summary tables are consistent")
        else:
            print("Field summary tables are inconsistent")
    return test_rslt


# define command line interface!
parser = argparse.ArgumentParser(
    description = ("Measures and compares field summary statistics from "
                   "Enzo-E results.")
)
subparsers = parser.add_subparsers(help='sub-command help', dest = 'command')

def _add_target_path_arg(parser):
    parser.add_argument(
        'target_path', action = 'store',
        help = ("The path to the block_list file of the simulation from which "
                "the field properties should be computed.")
    )

measure_parser = subparsers.add_parser(
    'measure', help = 'Measures the field summary information'
)
_add_target_path_arg(measure_parser)
measure_parser.add_argument(
    '-f','--fields', action = StoreCommaListAction, required = True,
    help = ("Specify the list of fields to use to compute the error norm. The "
            "field names should be separated by commas and have no spaces.")
)
measure_parser.add_argument(
    '-o','--output', required = True,
    help = ("Path to file where the table should be written. Passing a hyphen "
            "indicates that it should be written to the terminal")
)
measure_parser.set_defaults(func = _main_measure)



cmp_parser = subparsers.add_parser(
    'cmp', help = 'measure the field summary and compares against a reference'
)
_add_target_path_arg(cmp_parser)
cmp_parser.add_argument(
    '--ref', required = True, action = 'store',
    help = "Path to the file of reference field summary information"
)
cmp_parser.add_argument(
    '--report', default = None,
    help = ('Path to file where a properly formatted test report describing '
            'the sucess of the comparison should be written.')
)
cmp_parser.add_argument(
    '--rtol', required = True, action = 'store', type = float,
    help = "Relative tolerance"
)
cmp_parser.add_argument(
    '--atol', required = True, action = 'store', type = float,
    help = "Absolute tolerance"
)

cmp_parser.set_defaults(func = _main_cmp)

if __name__ == '__main__':
    args = parser.parse_args()
    main_func = args.func
    main_func(args)
    
    """
    #print(colnames)
    old_field_table = measure_field_summary(
        './GRACKLE_TEST_000/GRACKLE_TEST_000.block_list', field_names
    )

    print("Constructing Field Summary Table:")
    field_table = measure_field_summary(
        './GRACKLE_TEST_010/GRACKLE_TEST_010.block_list', field_names
    )
    print(old_field_table == field_table)
    #print(field_table)
    #print(field_table['min','max'])

    print("Comparing Field Summary Table:")

    with create_test_report('my_report.unit', clobber = True) as tr:
        test_rslt = test_equivalent_field_tables(old_field_table,field_table,
                                                 atol = 1e-15, rtol = 1e-14,
                                                 test_report = tr)
    """
