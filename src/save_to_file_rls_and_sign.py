__author__="huziy"
__date__ ="$Aug 29, 2011 12:46:18 PM$"



def save_to_file_rls_and_sign(current_id, return_period,
                              return_levels_current, return_levels_future,
                              rl_stdevs_current, rl_stdevs_future,
                              condition, current_highs, future_highs):

    outPath = '{0}_{1}_rl_and_sign.txt'.format(current_id, return_period)

    f = open(outPath, 'w')
    header = 'current, future \n'
    footer0 = 'return levels (current, future): {0}, {1} \n'
    footer1 = 'ret. lev stds (current, future): {0}, {1} \n'
    footer = 'significant: {0}\n\n'
    line_format = '{0},{1}\n'

    print(len(current_highs))
    npos = len(current_highs)
    for pos in range(npos):
#        f.write(header)

        the_highs = []
        the_highs.extend(current_highs[pos])
        the_highs.extend(future_highs[pos])


        the_highs = ['%.2f' % x for x in the_highs]
        f.write(';'.join(the_highs) + '\n')

#        for current, future in zip(current_highs[pos], future_highs[pos]):
#            f.write(line_format.format(current, future))
#        f.write(footer0.format(return_levels_current[pos], return_levels_future[pos]))
#        f.write(footer1.format(rl_stdevs_current[pos], rl_stdevs_future[pos]))
#        f.write(footer.format(condition[pos]))
#        f.write(20 * '-' + '\n')
    f.close()
    pass


if __name__ == "__main__":
    print("Hello World")
