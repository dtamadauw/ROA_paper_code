#TITLE: MATLAB code for "Motion-robust gradient-echo T2 mapping in the liver using readout alternation"
#Matlab & Python scripts to generate table S1 in the paper
#Author: Daiki Tamada
#Affiliation: Department of Radiology, University of Wisconsin-Madison
#Date: 9/21/2026
#Email: dtamada@wisc.edu

#By downloading, installing, or otherwise accessing or using the Software , you ("Recipient") agree to receive and use the above-
#identified SOFTWARE subject to the following terms, obligations and restrictions. If you do not agree to all of the following terms, 
#obligations and restrictions you are not permitted to download, install,
#execute, access, or use the SOFTWARE:

#1.	Originators of the SOFTWARE.  Provider is willing to license its rights in the SOFTWARE ("Provider's Rights") to academic researchers to use free of charge solely for academic, non-commercial research purposes subject to the terms and conditions outlined herein. The SOFTWARE was created at the University of Wisconsin ("UW") by Dakai Tamada. Please note Provider's Rights may include, but are not limited to, certain patents or patent applications owned by the Wisconsin Alumni Research Foundation ("WARF"). 
#2.	Limited License.  Provider hereby grants to Recipient a non-commercial, non-transferable, royalty-free, non-exclusive license, without the right to sublicense, under Provider's Rights to  download, install, access, execute and use the SOFTWARE solely for academic, non-commercial research purposes. SOFTWARE may not be used, directly or indirectly, to perform services for a fee or for the production or manufacture of products for sale to third parties. The foregoing license does not include any license to third party intellectual property that may be contained in the SOFTWARE; obtaining a license to such rights is Recipient's responsibility. 
#3.	Restrictions on SOFTWARE use and distribution.  Recipient shall not take, authorize, or permit any of the following actions with the SOFTWARE: (1) modify, translate or otherwise create any derivative works; or (2) publicly display (e.g., Internet) or publicly perform (e.g., present at a press conference); or (3) sell, lease, rent or lend; or (4) use it for any commercial purposes whatsoever. Recipient must fully reproduce and not obscure, alter or remove any of the Provider's proprietary notices that appear on the SOFTWARE, including copyright notices or additional license terms included with any the third party software contained in the SOFTWARE. Recipient may not provide any third party with access to the SOFTWARE or use the SOFTWARE on a timeshare or service bureau basis. Recipient represents that it is compliance with all applicable export control provisions and is not prohibited from receiving the SOFTWARE. 
#4.	Reservation of rights.  Provider retains all rights and title in the SOFTWARE, including without limitation all intellectual property rights (e.g., patent, copyright and trade secret rights) that may now or in the future exist in the SOFTWARE, regardless of form or medium. Provider retains ownership and all of Its rights in the SOFTWARE, including all of its intellectual property rights (e.g., patent, copyright and trade secret rights) that may now or in the future cover the SOFTWARE or any uses of the SOFTWARE, regardless of form or medium; title remains with Provider and the SOFTWARE is merely being loaned to Recipient for the specific purposes and under the specific restrictions stated herein. Nothing in this Agreement grants Recipient any additional rights to the SOFTWARE, any right to obtain any updates or new releases of the SOFTWARE, any commercial license for the SOFTWARE, or any other intellectual property owned or licensed by Provider. Provider has no obligation to provide any support, updates, or bug fixes.
#5.	Disclaimer of Warranty. PROVIDER IS PROVIDING THE SOFTWARE TO RECIPIENT ON AN "AS IS" BASIS. PROVIDER MAKES NO REPRESENTATIONS OR WARRANTIES CONCERNING THE SOFTWARE OR ANY OUTCOME THAT MAY BE OBTAINED BY USING THE SOFTWARE, AND EXPRESSLY DISCLAIMS ALL SUCH WARRANTIES, INCLUDING WITHOUT LIMITATION ANY EXPRESS OR IMPLIED WARRANTY OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, AND NON-INFRINGEMENT OF INTELLECTUAL PROPERTY RIGHTS. PROVIDER MAKES NO REMEDY THAT THE SOFTWARE WILL OPERATE ERROR FREE OR UNINTERRUPTED.
#6.	Limitation of Liability; Indemnity.  TO THE FULLEST EXTENT PERMITTED BY LAW, IN NO EVENT SHALL PROVIDER BE LIABLE TO RECIPIENT FOR ANY LOST PROFITS OR ANY DIRECT, INDIRECT, EXEMPLARY, PUNITIVE, INCIDENTAL OR CONSEQUENTIAL DAMAGES ARISING FROM THE SOFTWARE OR ITS USE. FURTHERMORE, IN NO EVENT WILL PROVIDER'S LIABILITY TO RECIPIENT EXCEED $100. PROVIDER HAS NO LIABILITY FOR ANY DECISION, ACT OR OMISSION MADE BY RECIPIENT AS A RESULT OF USE OF THE SOFTWARE. To the extent permitted by applicable law, Recipient agrees to indemnify, defend and hold harmless Provider, UW, and the SOFTWARE authors against all claims and expenses, including legal expenses and reasonable attorneys fees, arising from Recipient's use of the SOFTWARE.
#7.	No use of names/trademarks.  Recipient shall not use Provider's name, or the name of any author of the SOFTWARE or that of UW, in any manner without the prior written approval of the entity or person whose name is being used.
#8.	Termination.  Without prejudice to any other rights, Provider may terminate this Agreement if Recipient fails to comply with the terms of this Agreement for any reason. Upon termination for any reason, Recipient must immediately destroy all copies of the SOFTWARE in Recipient's possession, custody, or control.	

from __future__ import annotations
import argparse, re
from pathlib import Path
import numpy as np
import pandas as pd

RF_PATTERN = re.compile(r"RF(\d+):(F\+|F-|Z)")
M1_SCALE = 1e3


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument('--conventional', default='./steady_initialized_pathways_conventional.csv')
    p.add_argument('--roa', default='./steady_initialized_pathways_roa.csv')
    p.add_argument('--output-dir', default='./')
    p.add_argument('--n-dominant', type=int, default=12)
    p.add_argument('--cumulative-power', type=float, default=0.90)
    p.add_argument('--allow-conventional-sl', action='store_true',
                   help='Retain conventional M1_SL from input. By default it is forced to zero.')
    return p.parse_args()


def rf_states(s):
    return [(int(tr), st) for tr, st in RF_PATTERN.findall(str(s))]


def pathway_key(s):
    x = rf_states(s)
    return '|'.join(f'RF{tr}:{st}' for tr, st in x) if x else str(s)


def compact_pattern(s):
    x = rf_states(s)
    return ' → '.join(f'{st}_{tr}' for tr, st in x if st != 'Z') if x else str(s)


def simple_description(s):
    x = rf_states(s)
    trs = [tr for tr, st in x if st != 'Z']
    if not trs:
        return 'No transverse signal during traced window'
    if len(trs) == 1:
        return f'Rephased signal appearing during TR{trs[0]}'
    z_between = any(st == 'Z' for _, st in x[x.index((trs[0], next(st for tr,st in x if tr==trs[0]))):])
    # Reader-friendly description only; exact EPG pattern remains in full table.
    contiguous = all(b == a + 1 for a, b in zip(trs[:-1], trs[1:]))
    if contiguous:
        return f'Signal remains transverse during TR{trs[0]}-{trs[-1]}'
    return 'Signal rephases after evolution over multiple TR intervals'


def weighted_rms(v, w):
    v = np.asarray(v, float); w = np.asarray(w, float)
    return np.nan if w.sum() == 0 else float(np.sqrt(np.sum(w*v*v)/np.sum(w)))


def aggregate(df):
    keys = ['PathwayKey','CompactPattern','Description','EchoType']
    rows = []
    for k, g in df.groupby(keys, dropna=False, sort=False):
        amp = g.AmplitudeReal.to_numpy() + 1j*g.AmplitudeImag.to_numpy()
        a = amp.sum(); w = np.abs(amp)**2
        if w.sum() == 0: w = np.ones(len(g))
        def signed(col):
            vals = g[col].to_numpy(float)
            sg = np.sign(np.sum(w*vals)) or 1.0
            return sg*weighted_rms(vals,w)
        rows.append(dict(PathwayKey=k[0], CompactPattern=k[1], Description=k[2], EchoType=k[3],
                         Pathway=g.Pathway.iloc[0], NParentStates=len(g),
                         AmplitudeReal=a.real, AmplitudeImag=a.imag,
                         Magnitude=abs(a), Power=abs(a)**2,
                         M1_RO=signed('M1_RO'), M1_SL=signed('M1_SL'),
                         M1_RSS=weighted_rms(g.M1_RSS,w)))
    out = pd.DataFrame(rows)
    out['RelativeMagnitude_pct'] = 100*out.Magnitude/out.Magnitude.sum()
    out['RelativePower_pct'] = 100*out.Power/out.Power.sum()
    for c in ['M1_RO','M1_SL','M1_RSS']:
        out[c+'_scaled'] = M1_SCALE*out[c]
    return out


def load(path, acquisition, force_zero_sl):
    df = pd.read_csv(path)
    req = {'Pathway','EchoType','M1_RO','M1_SL','M1_RSS','AmplitudeReal','AmplitudeImag'}
    miss = req - set(df.columns)
    if miss: raise ValueError(f'{path} missing columns: {sorted(miss)}')
    if acquisition == 'Conventional' and force_zero_sl:
        df['M1_SL'] = 0.0
        df['M1_RSS'] = np.abs(df['M1_RO'])
    df['PathwayKey'] = df.Pathway.map(pathway_key)
    df['CompactPattern'] = df.Pathway.map(compact_pattern)
    df['Description'] = df.Pathway.map(simple_description)
    return aggregate(df)


def dominant(df, n, cumulative):
    x = df.sort_values('RelativePower_pct', ascending=False).copy()
    x['CumulativePower_pct'] = x.RelativePower_pct.cumsum()
    n_cum = np.searchsorted(x.CumulativePower_pct.to_numpy(), 100*cumulative, side='left') + 1
    x = x.head(max(n, n_cum)).copy()
    x.insert(0, 'Rank', np.arange(1, len(x)+1))
    return x


def full_table(x):
    return pd.DataFrame({
        'Rank':x.Rank,
        'Compact EPG pattern':x.CompactPattern,
        'Reader-friendly description':x.Description,
        'Steady-state parent states':x.NParentStates,
        'Relative magnitude (%)':x.RelativeMagnitude_pct,
        'Relative pathway power (%)':x.RelativePower_pct,
        'Effective M1,RO (10^-3 T s^2/m)':x.M1_RO_scaled,
        'Effective M1,SL (10^-3 T s^2/m)':x.M1_SL_scaled,
        'Effective M1,RSS (10^-3 T s^2/m)':x.M1_RSS_scaled,
        'Cumulative pathway power (%)':x.CumulativePower_pct,
    }).round(3)


def reader_table(x):
    return pd.DataFrame({
        'Rank':x.Rank,
        'Signal component':[f'Component {i}' for i in x.Rank],
        'Description':x.Description,
        'Relative pathway power (%)':x.RelativePower_pct,
        'Effective readout first moment (10^-3 T s^2/m)':x.M1_RO_scaled,
        'Effective slice first moment (10^-3 T s^2/m)':x.M1_SL_scaled,
        'Effective total first moment (10^-3 T s^2/m)':x.M1_RSS_scaled,
        'Cumulative pathway power (%)':x.CumulativePower_pct,
    }).round(3)


def ensemble(df, name):
    w = df.Power
    return {'Acquisition':name,
            'Number of unique traced histories':len(df),
            'Number of steady-state parent states':int(df.NParentStates.sum()),
            'M1,RO RMS power-weighted (10^-3 T s^2/m)':M1_SCALE*weighted_rms(df.M1_RO,w),
            'M1,SL RMS power-weighted (10^-3 T s^2/m)':M1_SCALE*weighted_rms(df.M1_SL,w),
            'M1,RSS RMS power-weighted (10^-3 T s^2/m)':M1_SCALE*weighted_rms(df.M1_RSS,w)}


def main():
    a = parse_args(); out = Path(a.output_dir); out.mkdir(parents=True, exist_ok=True)
    force_zero = not a.allow_conventional_sl
    c = load(a.conventional,'Conventional',force_zero)
    r = load(a.roa,'ROA',force_zero)
    dc = dominant(c,a.n_dominant,a.cumulative_power)
    dr = dominant(r,a.n_dominant,a.cumulative_power)

    reader_table(dc).to_csv(out/'table_s1_dominant_pathways_conventional.csv',index=False)
    reader_table(dr).to_csv(out/'table_s1_dominant_pathways_roa.csv',index=False)


if __name__ == '__main__':
    main()
