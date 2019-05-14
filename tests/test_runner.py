import sys
import json
import hashlib
from subprocess import Popen, PIPE


def execute_test(name, test):
    print("Running test %s" % name)
    command = test.get('command')
    exit_code = test.get('exit_code')
    print("Command: %s" % command)
    proc = Popen(command, stdout=PIPE, universal_newlines=True)
    for line in proc.stdout:
        print(line, end='')
    proc.stdout.close()
    exit_code = proc.wait()
    err_msg = ""
    passed = True
    for check in test.get('checks'):
        if check['type'] == 'checkExitCode':
            successful, msg = check_exit_code(exit_code, check['expected'])
        elif check['type'] == 'checkFile':
            successful, msg = check_file(check['filename'], check['checksum'])
        elif check['type'] == 'checkNumberOfLines':
            successful, msg = check_number_of_lines(check['filename'], check['num_lines']) 
        if not successful:
                passed = False
                err_msg += "%s " % msg
    return passed, err_msg


def check_number_of_lines(filename, num):
    num_lines = sum(1 for line in open(filename))
    if num_lines == num:
        return True, ""
    else:
        return False, "Expected number of lines %s, got %s" % (num, num_lines)


def check_file_checksum(filename, checksum):
    with open(filename,"rb") as f:
        bytes = f.read() # read entire file as bytes
        readable_hash = hashlib.sha256(bytes).hexdigest();
    if readable_hash == checksum:
        return True, ""
    else:
        return False, "Expected checksum %s, got %s" % (checksum, readable_hash)


def check_exit_code(got, expected):
    if got == expected:
        return True, ""
    else:
        return False, "Expected exit_code %s, got %s" % (expected, got)


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
        passed, msg = execute_test(k, v)
        if passed:
            report.append("PASSED: %s" % k)
        else:
            report.append("FAILED: %s ERROR: %s" % (k, msg))
            successful = False
    print_report(report)
    if not successful:
        exit(1)
