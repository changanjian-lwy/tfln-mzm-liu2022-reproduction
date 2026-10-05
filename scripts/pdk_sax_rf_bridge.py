"""Tool-learning demo (not registered evidence): the A04 L1 electrode as a SAX S-parameter netlist.

Builds source-reference -> uniform line (Z0(f), gamma(f), L) -> load as a SAX circuit and compares
the input S11 with the project's closed form tfln_mzm.termination.s11. Then inserts a lumped pad
capacitance to show how an L3 (SC-02) network would be described; that part uses an illustrative
value only and is not a model input.
Needs the optional PDK environment (Python 3.12, requirements-pdk-lock.txt); see docs/PIC_DESIGN_FLOW_zh.md.
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import sax

REPO=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(REPO/'src'),str(REPO/'scripts')]
from run_b01 import inputs
from tfln_mzm.microwave import propagation
from tfln_mzm.termination import s11 as s11_closed

ZR=50.0


def line(f_GHz=1.0,z0=50.0,gamma=0j,length_m=1e-3):
    """Uniform TEM line referenced to ZR on both ports."""
    g=gamma*length_m;sh,ch=np.sinh(g),np.cosh(g)
    den=2*z0*ZR*ch+(z0**2+ZR**2)*sh
    s11=(z0**2-ZR**2)*sh/den;s21=2*z0*ZR/den
    return sax.reciprocal({('in','in'):s11,('out','out'):s11,('in','out'):s21})


def load(zl=50.0):
    return {('p','p'):(zl-ZR)/(zl+ZR)}


def shunt_c(f_GHz=1.0,c_f=0.0):
    """Shunt capacitor between two ports (illustrative pad parasitic)."""
    y=2j*np.pi*f_GHz*1e9*c_f*ZR
    s11=-y/(2+y);s21=2/(2+y)
    return sax.reciprocal({('in','in'):s11,('out','out'):s11,('in','out'):s21})


def main():
    d=pd.read_csv(REPO/'data/digitized/liu2022/curves.csv',dtype={'series':str})
    f,a,z0=inputs(d,np.arange(.025,200+1e-9,.025),False)
    gam=propagation(f*1e9,2.25,a)
    electrode,_=sax.circuit(netlist={'instances':{'tl':'line','zl':'load'},'connections':{'tl,out':'zl,p'},'ports':{'rf_in':'tl,in'}},
                            models={'line':line,'load':load})
    out={}
    for zl in (20.0,40.0,50.0,80.0):
        s_sax=np.asarray(electrode(tl={'z0':z0,'gamma':gam,'length_m':.006},zl={'zl':zl})['rf_in','rf_in'])
        s_ref=s11_closed(gam,.006,z0,zl)
        out[zl]=float(np.abs(s_sax-s_ref).max())
    print('max |S11_SAX - S11_closed| over',len(f),'supported points:',out)
    padded,_=sax.circuit(netlist={'instances':{'pad':'shunt_c','tl':'line','zl':'load'},
                                  'connections':{'pad,out':'tl,in','tl,out':'zl,p'},'ports':{'rf_in':'pad,in'}},
                         models={'shunt_c':shunt_c,'line':line,'load':load})
    s_pad=np.asarray(padded(pad={'f_GHz':f,'c_f':20e-15},tl={'z0':z0,'gamma':gam,'length_m':.006},zl={'zl':40.0})['rf_in','rf_in'])
    s0=s11_closed(gam,.006,z0,40.0);k=np.argmin(abs(f-100))
    print('illustrative 20 fF pad, ZL=40: max S11 0-200 GHz %.2f dB vs %.2f dB without; at %.1f GHz %.2f vs %.2f dB'%(
        20*np.log10(abs(s_pad).max()),20*np.log10(abs(s0).max()),f[k],20*np.log10(abs(s_pad[k])),20*np.log10(abs(s0[k]))))


if __name__=='__main__':
    main()
