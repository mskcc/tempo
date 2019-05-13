import sys
import json
from subprocess import Popen, PIPE


def execute_test(name, test):
    print("Running test %s" % name)
    command = test.get('command')
    print("Command: %s" % command)
    proc = Popen(command, stdout=PIPE, universal_newlines=True)
    for line in proc.stdout:
        print(line, end='')
    proc.stdout.close()
    exit_code = proc.wait()
    return exit_code


def print_report(status):
    REPORT = "\nTEST REPORT:\n------------------------------------------"
    for test in status:
        REPORT += "\n%s" % test
    print(REPORT)


if __name__ == '__main__':
    test_file = sys.argv[1]
    successful = True
    with open(test_file, 'r') as f:
        tests = json.load(f)
    report = []
    for k, v in tests.items():
        return_code = execute_test(k, v)
        if return_code == 0:
            report.append("PASSED: %s" % k)
        else:
            report.append("FAILED: %s" % k)
            successful = False
    print_report(report)
    if not successful:
        exit(1)
