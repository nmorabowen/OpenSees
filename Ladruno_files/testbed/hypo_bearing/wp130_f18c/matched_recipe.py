"""WP-130 round 1: the recommended recipe's committed q-s vs IntScheme 1, matched s/B."""
import analyse as A
import matched as M

ARMS = ['s1_T2', 'rec_s1_T2', 'rec_s2']


def rel(q, r):
    return '--' if (q is None or r is None) else f'{abs(q - r) / r * 100:.2f} %'


print('| s/B | ' + ' | '.join(ARMS) + ' | rec_s2 vs s1_T2 | rec_s2 vs rec_s1_T2 |')
print('|---|' + '---|' * (len(ARMS) + 2))
for sb in (0.001, 0.002, 0.004, 0.006):
    qs = [M.q_at(A.curve(a), sb) for a in ARMS]
    print(f'| {sb} | ' + ' | '.join('--' if q is None else f'{q:.2f}' for q in qs)
          + f' | {rel(qs[2], qs[0])} | {rel(qs[2], qs[1])} |')
