import matplotlib.pyplot as plt
import numpy as np
import math #needed for definition of pi
import pandas as pd
import csv

damage_coarse = pd.read_csv('iso_calibration_continuum_coarse_Reaction.txt', sep = "   ", names = ['displacement','reaction_force'])
damage_coarse.displacement = abs(damage_coarse.displacement)
damage_coarse.reaction_force = abs(damage_coarse.reaction_force)

damage_fine = pd.read_csv('iso_calibration_continuum_fine_Reaction.txt', sep = "   ", names = ['displacement','reaction_force'])
damage_fine.displacement = abs(damage_fine.displacement)
damage_fine.reaction_force = abs(damage_fine.reaction_force)

damage_czm = pd.read_csv('iso_calibration_CZM_Reaction.txt', sep = "   ", names = ['displacement','reaction_force'])
damage_czm.displacement = abs(damage_czm.displacement)
damage_czm.reaction_force = abs(damage_czm.reaction_force)

plt.figure(1)
plt.plot(damage_coarse.displacement,damage_coarse.reaction_force, label = "damage_coarse")
plt.plot(damage_fine.displacement,damage_fine.reaction_force, label = "damage_fine")
plt.plot(damage_czm.displacement,damage_czm.reaction_force, label = "CZM")

# plt.plot(year, prod_pp,label = "production well")
plt.xlabel("Displacement, mm")
plt.ylabel("Reaction force, N")
# plt.title('evolution')
plt.legend()

# # plt.figure(2)
# plt.subplot(212)
# plt.plot(year,inj_T, label = "injection well")
# plt.plot(year, prod_T,label = "production well")
# plt.xlabel("year")
# plt.ylabel("Temperature")
# # plt.title('evolution')
# plt.legend()
plt.show()


# with open('/Users/jinw-mac/projects/falcon/Applied_Energy/Cranfield_HT_season_CSV.csv', mode='r') as csv_file:
#     csv_reader = csv.reader(csv_file, delimiter=',')
#     line_count = 0

# x = np.arange(0, math.pi*2, 0.05)
# y = np.tan(x)
# 

