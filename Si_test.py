from __future__ import division
from BandDiagrams import *

from TunnelCalcs import *

# Generate a band diagram of straight lines for a silicon PIN junction
figure()
si_flatpin=generate_flatpin_bandstructure("Si",1e-5,"B",1e20,2e-6,1e-5,"P",1e19,linspace(-10,-.5,10))
s=si_flatpin[0]
plot(s["x"]*1e7,s["Ev"])
plot(s["x"]*1e7,s["Ec"])
plot(s["x"]*1e7,s["EB"],'--')
ylabel('Energy [eV]')
xlabel('Position [nm]')

# Generate a depletion approximation band diagram for a silicon PIN junction
figure()
si_pin=generate_pin_bandstructure("Si",1e-5,"B",1e20,2e-6,1e-5,"P",1e19,linspace(-10,-.5,10))
s=si_pin[0]
plot(s["x"]*1e7,s["Ev"])
plot(s["x"]*1e7,s["Ec"])
plot(s["x"]*1e7,s["EB"],'--')
ylabel('Energy [eV]')
xlabel('Position [nm]')

bsf=si_flatpin
jsf=[current(b,"orig") for b in bsf]
jsfd=[current(b,"v2") for b in bsf]
jsfm=[current(b,"maxfield") for b in bsf]
jsfa=[current(b,"meanfield") for b in bsf]
bsc=si_pin
jsc=[current(b,"orig") for b in bsc]
jscd=[current(b,"v2") for b in bsc]
jscm=[current(b,"maxfield") for b in bsc]
jsca=[current(b,"meanfield") for b in bsc]




figure(figsize=(8,6))
plot([max_field(b)/1e6 for b in bsf],
     [j*1e3/1e12/b["Vr"] for j,b in zip(jsf,bsf)],'b',label='Integral')
plot([max_field(b)/1e6 for b in bsf],
     [j*1e3/1e12/b["Vr"] for j,b in zip(jsfm,bsf)],'r',label='Analytic Max')
plot([max_field(b)/1e6 for b in bsf],
     [j*1e3/1e12/b["Vr"] for j,b in zip(jsfa,bsf)],'g',label='Analytic Mean')


xlim([0,4])
xlabel('Max Field [MV/cm]')
yscale('log')
gca().set_xticks([0,1,2,3,4])
gca().set_yticks([1e-12,1e-10,1e-8,1e-6,1e-4,1e-2,1,1e2])
ylim([1e-12,1e3])
ylabel('J/V [mA/um$^2$/V]')
legend(loc='upper left')
title('Linear PIN')

figure(figsize=(8,6))
plot([max_field(b)/1e6 for b in bsc],
     [j*1e3/1e12/b["Vr"] for j,b in zip(jsc,bsc)],'b.',label='Integral')
plot([max_field(b)/1e6 for b in bsc],
     [j*1e3/1e12/b["Vr"] for j,b in zip(jscm,bsc)],'r.',label='Analytic Max')
plot([max_field(b)/1e6 for b in bsc],
     [j*1e3/1e12/b["Vr"] for j,b in zip(jsca,bsc)],'g.',label='Analytic Mean')

xlim([0,4])
xlabel('Max Field [MV/cm]')
yscale('log')
gca().set_xticks([0,1,2,3,4])
gca().set_yticks([1e-12,1e-10,1e-8,1e-6,1e-4,1e-2,1,1e2])
ylim([1e-12,1e3])
ylabel('J/V [mA/um$^2$/V]')
legend(loc='upper left')
title('Quadratic PIN')

figure(figsize=(8,6))
plot([max_field(b)/1e6 for b in bsf],
     [j*1e3/1e12/b["Vr"]for j,b in zip(jsf,bsf)],'ws',label='Linear',markersize=10)
plot([max_field(b)/1e6 for b in bsc],
     [j*1e3/1e12/b["Vr"] for j,b in zip(jsc,bsc)],'b.',label='Integral',markersize=10)
plot([max_field(b)/1e6 for b in bsc],
     [j*1e3/1e12/b["Vr"] for j,b in zip(jscm,bsc)],'r.',label='Analytic Max',markersize=10)
plot([max_field(b)/1e6 for b in bsc],
     [j*1e3/1e12/b["Vr"] for j,b in zip(jsca,bsc)],'g.',label='Analytic Mean',markersize=10)

xlim([0,4])
xlabel('Max Field [MV/cm]')
yscale('log')
gca().set_xticks([0,1,2,3,4])
gca().set_yticks([1e-12,1e-10,1e-8,1e-6,1e-4,1e-2,1,1e2])
ylim([1e-12,1e3])
ylabel('J/V [mA/um$^2$/V]')
title('Linear vs Quadratic PIN')
legend(loc='upper left')


figure(figsize=(8,6))
plot([b["Vr"] for b in bsf],
     [j*1e3/1e12/b["Vr"] for j,b in zip(jsf,bsf)],linestyle='None',color='w',marker='s',markersize=10,label='Linear')
plot([b["Vr"] for b in bsc],
     [j*1e3/1e12/b["Vr"] for j,b in zip(jsc,bsc)],'b.',label='Integral',markersize=10)
plot([b["Vr"] for b in bsc],
     [j*1e3/1e12/b["Vr"] for j,b in zip(jscm,bsc)],'r.',label='Analytic Max',markersize=10)
plot([b["Vr"] for b in bsc],
     [j*1e3/1e12/b["Vr"] for j,b in zip(jsca,bsc)],'g.',label='Analytic Mean',markersize=10)

xlabel('Reverse bias [V]')
yscale('log')
gca().set_yticks([1e-12,1e-10,1e-8,1e-6,1e-4,1e-2,1,1e2])
ylim([1e-12,1e3])
ylabel('J/V [mA/um$^2$/V]')
title('Linear vs Quadratic PIN')
legend(loc='upper left')