from datetime import datetime

def print_header(terminal_width=80,
                 author='TJR, 2021',
                 title='COMAP Aperture Photometry Program',
                 version=None):
    '''
    Print a neatly formatted header to the terminal comprising of author, date &
    time, title and version information.
    '''

    now = datetime.now()
    d_string = now.strftime("%Y-%m-%d")
    dt_string = now.strftime("%Y-%m-%d_%H:%M:%S")

    if version != None:
        title = "{} (v{})".format(title,version)

    empty_2 = terminal_width-(len(dt_string) + len(author))
    empty_1 = terminal_width-(len(title))

    print('='*terminal_width)
    print(' '*(empty_1//2) + title)
    print('='*terminal_width)
    print(dt_string + ' '*empty_2 + author + '\n')


def parse_file(param_file_name):
    '''
    Parse input files into a dictionary for use in general purpose programs with
    specific modular pipelines. File should be formatted as below:

    [category1]
    {cat1_instance 1}
    cat1_ins1_var1  | 0
    cat1_ins1_var2  | 1
    ...
    {instance n}
    cat1_insn_var1  | 3
    cat1_insn_var2  | 4

    [category2]
    ...
    '''

    param_file = open(param_file_name, 'r')

    params = {}
    obj, obj2, var, val = None, None, None, None

    for line in param_file:

        if (line[0]!='#') and (line.strip() != ''):
            line = line[:-1].strip()

            if line[0] == '[':
                obj = line[1:].split(']')[0]
                params[obj] = {}
                obj2, var, val = None, None, None

            elif line[0] == '{':
                obj2 = line[1:].split('}')[0]
                params[obj][obj2] = {}
                var, val = None, None

            else:
                var, val = line.split('|')
                var = var.strip()

                # split val into list of values (or keep as one if "...")
                if val.strip()[0] == '"':
                    val = [val.strip()[1:-1]]
                else:
                    val = val.split()

                # get and apply datatype (int, float, frac, string)
                for i in range(len(val)):
                    alpha_flag, div_flag, dp_flag = False, False, False
                    for j in val[i]:
                        if j.isalpha():
                            alpha_flag = True
                        if j == "/":
                            div_flag = True
                        if j == ".":
                            dp_flag = True
                    if alpha_flag == True and div_flag == False and dp_flag == False:
                        if val[i].strip() == 'False':
                            val[i] = False
                        if val[i].strip() == 'True':
                            val[i] = True

                    if alpha_flag == False:
                        if div_flag == True:
                            temp = val[i].split('/')
                            val[i] = float(temp[0]) / float(temp[1])
                        else:
                            if dp_flag == True:
                                val[i] = float(val[i])
                            else:
                                val[i] = int(val[i])
                if len(val) == 1:
                    val = val[0]

                # write out to dictionary
                if obj2 != None:
                    params[obj][obj2][var] = val
                else:
                    params[obj][var] = val
    return params
