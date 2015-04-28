"""This module can be used to generate enumerated progress messages."""

import re

prog_ct = 0
sub_ct = 0

def progress(msg, done=True):
    """Prints an enumerated progress message

    Keyword arguments:
    msg  -- string to be printed
    done -- done message will be used
    """

    global prog_ct

    print('[' + str(prog_ct) + ']', msg + '...')
    if not done:
        prog_ct += 1

def subprogress_init_counter():
    """Initialiye counter for subprogress"""

    global sub_ct
    sub_ct = 0

def subprogress(num, msg='', counter=True, fn_filter=(lambda x: x)):
    """Print subprogress message with progress percentage.

    You can use %% inside the message string to insert the current counter value
    into the message.

    Keyword arguments:
    num       -- number of all subprogress messages
    msg       -- string to be printed
    counter   -- should counter be increased
    fn_filter -- changes counter value if inserted into the message string
    """

    global prog_ct, sub_ct

    s = ''
    new_msg =''

    if counter:
        s = str(round(sub_ct / num * 100, 2)) + '%'

        m = re.split('%%', msg)

        for substr in m[:-1]:
            new_msg += substr
            new_msg += str(fn_filter(sub_ct))

        new_msg += m[-1]
    else:
        s = str(num)

    print('[' + str(prog_ct) + '][' + s + ']', new_msg)

    sub_ct += 1

def done(msg=''):
    """Print an enumerated done message and increase the counter

    msg -- can be used to comment what is finished here
    """

    global prog_ct

    print('[' + str(prog_ct) + ']', 'Done.')
    prog_ct += 1

def all_done():
    """Tell all that the program is finished"""

    global prog_ct

    print('[' + str(prog_ct) + ']', 'Finished all.')
