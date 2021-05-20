import contextlib
import os.path
import sys
import traceback

# Encapsulates test result reporting. The interface is probably not optimal,
# but at least it's a start.
#
# The interface is loosely modelled off of python's file API.

class TestReport:
    """
    Respresents the report for the current test(s).

    DO NOT CONSTRUCT THIS DIRECTLY. USE create_test_report INSTEAD.

    Note
    ----
    The underlying testing infrastructure is a little confusing. After looking
    through things, tests results are read from individual files. 

    The test file must include a line that says "UNIT TEST BEGIN". If it 
    doesn't include "END CELLO" and "UNIT TEST END", anywhere, the test is 
    counted as crashed (this count seems to happen in 2 different places). 
    This doesn't actually affect the end result of the test suite.

    Each test can report multiple test cases that are failures, are incomplete,
    or pass. The test suite (and CircleCi) only fails if at least one failure 
    is reported

    The underlying testing infrastructure requires that this is at least 
    somewhat eager about writing information to file. For example, if there is 
    a test failure, but some (related or unrelated) problem causes the program
    this to exit before writing that to file, the test infrastructure will 
    currently indicate success.
    """
    def __init__(self, f):
        self._f = f
        self._complete = False
        self.write('UNIT TEST BEGIN\n',flush = True)

    def write(self, msg, flush = False):
        self._f.write(msg)
        if flush:
            self._f.flush()

    def _test_case_status(self, prefix, msg = None):
        # it's EXTREMELY important that there is a space between the newline
        # and the prefix
        prefix = '\n ' + prefix
        if msg is not None:
            msg = ': '.join([prefix,msg])
        else:
            msg = prefix
        self.write(msg, flush = True)

    def fail(self, msg = None):
        # don't touch the prefix string unless it changes in build.sh
        self._test_case_status(prefix = "FAIL", msg = msg)

    def incomplete(self, msg = None):
        # If the reason you are invoking this method should prevent continuous
        # integration from passing, you should use `fail` instead of method

        # don't touch the prefix string unless it changes in build.sh
        self._test_case_status(prefix = "incomplete", msg = msg)

    def passing(self,msg = None):
        # don't touch the prefix string unless it changes in build.sh
        self._test_case_status(prefix = "pass", msg = msg)

    def complete(self):
        # At the time of writing, test results reported after complete will
        # still be counted. It's also currently OK for test reports to call
        # this method more than once. However, it's unclear whether that could
        # change in the future
        self._complete = True
        self.write('\nUNIT TEST END'
                   '\nEND CELLO\n',
                   flush = True)

    def is_complete(self):
        return self._complete

@contextlib.contextmanager
def create_test_report(test_file, clobber = True):
    """
    This is a context manager used to construct an instance of TestReport.

    If the TestReport has not already completed at the end of the `with` 
    statement, this will call it's `complete` method unless there is an 
    exception - in that case, `fail` will be called and the report will be left
    incomplete. If an exception but the test report is already `complete`,
    the exception will just be printed
    """

    if os.path.isfile(test_file):
        if clobber:
            print("Overwriting {}".format(test_file))
        else:
            raise ValueError(
                "The test report file, {}, already exists".format(test_file)
            )
    else:
        print("Creating {}".format(test_file))
    f = open(test_file,'w')
    test_report = TestReport(f)

    err = None
    err_str = None
    try:
        yield test_report
    except:
        exc_info = sys.exc_info()
        err = exc_info[1]
        err_str = ''.join(traceback.format_exception(*exc_info))
    finally:
        if not test_report.is_complete():
            if err_str is not None:
                test_report.fail(
                    "Unexpected python exception occured before the test "
                    "completed. Exception information is provided below:\n"
                )
                test_report.write(err_str)
                test_report.write('\n')
                # leave the report incomplete
            else:
                test_report.complete()
        elif err_str is not None:
            test_report.write(
                "Unexpected python exception occured after the test "
                "completed. Exception information is provided below:\n"
            )
            test_report.write(err_str)
            test_report.write('\n')

        f.close()

    # if an error occured, let's pass it through
    if err is not None:
        raise err
