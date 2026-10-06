#!/usr/bin/env python3
"""Write an NNLOJET runcard for one part from the ZEUS template
(nnlo21/validation/nnlojet_epLJJ_zeus2j.run): the placeholders SEED,
WARMUP/PRODUCTION become one sweep line, and the CHANNELS block gets the
part (LO, R, V as whole channels; luminosity parts such as RRa_3, RV_6 from
`NNLOJET -listlumi epLJJ`, with the region of RRa/RRb, as NNLOJET's own
workflow dokan does).

  nnlojet_card.py TEMPLATE PART warmup|production NCALL NITER SEED > job.run
"""
import re, subprocess, sys, os

NNLOJET = os.environ.get('NNLOJET', '/ptmp/mpp/akarlber/disorder-nnlo21/nnlojet/bin/NNLOJET')


def lumi(part):
    out = subprocess.run([NNLOJET, '-listlumi', 'epLJJ'], capture_output=True, text=True, check=True).stdout
    for l in out.splitlines():
        m = re.match(r'^\s*(\w+)\s+(.*! channel:.*)$', l)
        if m and m.group(1) == part:
            return m.group(2)
    sys.exit('no luminosity part ' + part)


def main():
    tmpl, part, mode, ncall, niter, seed = sys.argv[1:7]
    if part in ('LO', 'R', 'V', 'RR', 'RV', 'VV'):
        head, chans = 'CHANNELS', part
    else:
        m = re.match(r'^(RR)([ab])_\d+$', part)
        head = 'CHANNELS  region = %s' % m.group(2) if m else 'CHANNELS'
        chans = lumi(part)
    txt = open(tmpl).read()
    txt = re.sub(r'iseed = SEED', 'iseed = %s' % seed, txt)
    txt = re.sub(r'\n\s*warmup = WARMUP\s*\n', '\n  %s = %s[%s]\n' % (mode, ncall, niter), txt)
    txt = re.sub(r'\n\s*production = PRODUCTION\s*\n', '\n', txt)
    txt = re.sub(r'\nCHANNELS\s*\n\s*CHANNEL\s*\nEND_CHANNELS', '\n%s\n  %s\nEND_CHANNELS' % (head, chans), txt)
    assert 'SEED' not in txt.split('!')[0] and '\n  CHANNEL\n' not in txt
    sys.stdout.write(txt)


main()
