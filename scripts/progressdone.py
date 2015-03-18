import re

prog_ct = 0
sub_ct = 0

def progress(msg, done=True):
    global prog_ct

    print('[' + str(prog_ct) + ']', msg + '...')
    if not done:
            prog_ct += 1

def subprogress_init_counter():
    global sub_ct
    sub_ct = 0

def subprogress(num, msg='', counter=True):
    global prog_ct, sub_ct

    s = ''
    new_msg =''

    if counter:
        s = str(sub_ct / num * 100) + '%'

        m = re.split('%%', msg)

        for substr in m[:-1]:
            new_msg += substr
            new_msg += str(sub_ct)

        new_msg += m[-1]
    else:
        s = str(num)

    print('[' + str(prog_ct) + '][' + s + ']', new_msg)

    sub_ct += 1

def done(msg=''):
    global prog_ct

    print('[' + str(prog_ct) + ']', 'Done.')
    prog_ct += 1

def all_done():
    global prog_ct

    print('[' + str(prog_ct) + ']', 'Finished all.')
