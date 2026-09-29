"""WP-130: committed load-settlement at matched s/B, IntScheme-2 arms vs s1_T2."""
import analyse as A

ARMS = ['s1_T2', 's2t_T2', 's2t_ref0', 's2t_ref0_se_ls', 's2t_ref3_se_ls']


def q_at(rows, sb):
    prev = None
    for r in rows:
        x = float(r['s_over_B'])
        if x >= sb:
            if prev is None:
                return None
            x0, q0 = float(prev['s_over_B']), float(prev['q_foot_kPa'])
            return q0 + (sb - x0) / (x - x0) * (float(r['q_foot_kPa']) - q0)
        prev = r
    return None


print('| s/B | ' + ' | '.join(ARMS) + ' | max rel. diff vs s1_T2 |')
print('|---|' + '---|' * (len(ARMS) + 1))
for sb in (0.001, 0.002, 0.004, 0.006, 0.008):
    qs = [q_at(A.curve(a), sb) for a in ARMS]
    ref = qs[0]
    ds = [abs(q - ref) / ref for q in qs[1:] if q is not None and ref]
    d = max(ds) if ds else float('nan')
    print(f'| {sb} | ' + ' | '.join('--' if q is None else f'{q:.2f}' for q in qs)
          + f' | {d * 100:.2f} % |')
