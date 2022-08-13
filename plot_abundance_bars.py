import numpy as np
import matplotlib.pyplot as plt 

######################################################################## Creating stacked bar graph of GGchem abundances #################################################################################

# Zone 1 
# ~ abund_z1 = {'Fe': 10**-4.5049278777, 'Mg2SiO4': 10**-4.7609480571, 'SiO': 10**-4.948601267, 'Ni': 10**-5.8131482794, 'MgAl2O4': 10**-5.8510432007}        # 0.19076475972781792 AU
abund_z1 = {'Fe': 10**-4.5000000262, 'Mg2SiO4': 10**-4.8436701727, 'MgSiO3': 10**-5.0205869624, 'Ni': 10**-5.7800000844, 'NaAlSi3O8': 10**-5.8535825859}        # 0.39048869326270547 AU

abund_z1_total = 0
for val in abund_z1.values():
	abund_z1_total += val
abund_z1_percent = {key: value / abund_z1_total for key, value in abund_z1.items()}
print(abund_z1_percent)

# Zone 2
abund_z2 = {'Fe': 10**-4.5000000001, 'Mg2SiO4': 10**-4.8241760612, 'MgSiO3': 10**-5.074068744, 'NaAlSi3O8': 10**-5.7615469841, 'Ni': 10**-5.78}              #  1.0019278131870175 AU
abund_z2_total = 0
for val in abund_z2.values():
	abund_z2_total += val
abund_z2_percent = {key: value / abund_z2_total for key, value in abund_z2.items()}
print(abund_z2_percent)

# Zone 3
abund_z3 = {'Mg2SiO4': 10**-4.7580521304, 'FeS': 10**-4.9180726814, 'Fe2SiO4': 10**-5.0864970967, 'NaMg3AlSi3O12H2': 10**-5.8472306019, 'Fe3O4': 10**-5.9900858342}        # 10.127885287176355 AU
abund_z3_total = 0
for val in abund_z3.values():
	abund_z3_total += val
abund_z3_percent = {key: value / abund_z3_total for key, value in abund_z3.items()}
print(abund_z3_percent)

top1 = np.array([abund_z1_percent['Fe'], abund_z2_percent['Fe'], abund_z3_percent['FeS']])
top2 = np.array([abund_z1_percent['Mg2SiO4'], abund_z2_percent['Mg2SiO4'], abund_z3_percent['Mg2SiO4']])
top3 = np.array([abund_z1_percent['MgSiO3'], abund_z2_percent['MgSiO3'], abund_z3_percent['Fe2SiO4']])
top4 = np.array([abund_z1_percent['Ni'], abund_z2_percent['NaAlSi3O8'], abund_z3_percent['NaMg3AlSi3O12H2']])
top5 = np.array([abund_z1_percent['NaAlSi3O8'], abund_z2_percent['Ni'], abund_z3_percent['Fe3O4']])

# Plotting the bar chart
barWidth = 0.85
r = [0, 1, 2]
names = ['Zone 1', 'Zone 2', 'Zone 3']

fig, ax = plt.subplots(figsize=(8,15))
# ~ p1 = ax.bar(r, top1, color=['#b5ffb9', '#b5ffb9', '#a3acff'], edgecolor='white', width=barWidth)
# ~ p2 = ax.bar(r, top2, bottom = top1, color='#f9bc86', edgecolor='white', width=barWidth)
# ~ p3 = ax.bar(r, top3, bottom = top1+top2, color=['plum', 'plum', 'silver'], edgecolor='white', width=barWidth)
# ~ p4 = ax.bar(r, top4, bottom = top1+top2+top3, color=['paleturquoise', 'lightpink', 'lightcoral'], edgecolor='white', width=barWidth)
# ~ p5 = ax.bar(r, top5, bottom = top1+top2+top3+top4, color=['lightpink', 'paleturquoise', 'mediumpurple'], edgecolor='white', width=barWidth)

p1 = ax.bar(r, top2, color='#a3acff', edgecolor='white', width=barWidth)
p3 = ax.bar(r, top4, bottom = top2, color=['#b5ffb9', 'lightpink', 'mediumpurple'], edgecolor='white', width=barWidth)
p4 = ax.bar(r, top5, bottom = top2+top4, color=['lightpink', '#b5ffb9', 'lightcoral'], edgecolor='white', width=barWidth)
p4 = ax.bar(r, top3, bottom = top2+top4+top5, color=['#f9bc86', '#f9bc86', 'silver'], edgecolor='white', width=barWidth)
p5 = ax.bar(r, top1, bottom = top2+top4+top5+top3, color=['paleturquoise', 'paleturquoise', 'lemonchiffon'], edgecolor='white', width=barWidth)
ax.set_xticks(r, names)
# ~ ax.bar_label(p1, label_type='center')
 
# Custom x axis
ax.set_xticks(r, names)#, fontsize=14)
ax.set_xlabel("Radii")#, fontsize=16)
ax.set_ylabel("Mineral Abundances")#, fontsize=16)
fig.tight_layout()
plt.savefig("abundance_bars.png")
plt.show()
