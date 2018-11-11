import numpy as np
from math import *


'''
torques = M * acc + C * vel + N
'''
class dynamic():
    def __init__(self, dof, delta_t = 0.02):
        self.dof = dof
        self.M = np.zeros((self.dof,self.dof))
        self.C = np.zeros((self.dof,self.dof))
        self.N = np.zeros(self.dof)
        self.pos = np.zeros(self.dof)
        self.vel = np.zeros(self.dof)
        self.acc = np.zeros(self.dof)
        self.delta_t = delta_t


        self.lbx = 0
        self.lbz = 0
        self.l1x = 0.024
        self.l1z = 0.115
        self.l2x = 0.033
        self.l2z = 0
        self.l3x = 0
        self.l3z = 0.155
        self.l4x = 0
        self.l4z = 0.135
        self.l5x = 0
        self.l5z = 0.1136
        self.lfx = 0
        self.lfz = 0.05716

        self.l1m = 2.351
        self.l2m = 1.318
        self.l3m = 0.821
        self.l4m = 0.769
        self.l5m = 0.687
        self.mg = 0.219

        self.c1x = 0.01516
        self.c1y = 0.00359
        self.c1z = -0.03105
        self.c2x = -0.01903
        self.c2y = 0.0150
        self.c2z = 0.11397
        self.c3x = 0.00013
        self.c3y = 0.02022
        self.c3z = 0.10441
        self.c4x = 0.00015
        self.c4y = -0.02464
        self.c4z = 0.05353
        self.c5x = 0
        self.c5y = 0.0012
        self.c5z = 0.01648
        self.cgx = 0
        self.cgy = 0
        self.cgz = 0.0289

        self.I1x = 0.0029525
        self.I1y = 0.0060091
        self.I1z = 0.0058821
        self.I2x = 0.0031145
        self.I2y = 0.0005843
        self.I2z = 0.0031631
        self.I3x = 0.00172767
        self.I3y = 0.00041967
        self.I3z = 0.0018468
        self.I4x = 0.0006764
        self.I4y = 0.0010573
        self.I4z = 0.0006610
        self.I5x = 0.0001934
        self.I5y = 0.0001602
        self.I5z = 0.0000689
        self.Igx = 0.0002324
        self.Igy = 0.0003629
        self.Igz = 0.0002067


    def transition(self, torque):
        self.M_matrix()
        self.C_matrix()
        self.N_vector()
        self.acc = np.dot(np.linalg.pinv(self.M), (torque - self.C * self.vel - self.N))
        self.pos = self.pos + self.vel * self.delta_t + 0.5 * self.acc * self.delta_t **2
        self.vel = self.vel + self.acc * self.delta_t

        return self.pos, self.vel


    def M_matrix(self):
        pos = self.pos
        M = self.M
        lbx = 0;
        lbz = 0;
        l1x = 0.024;
        l1z = 0.115;
        l2x = 0.033;
        l2z = 0;
        l3x = 0;
        l3z = 0.155;
        l4x = 0;
        l4z = 0.135;
        l5x = 0;
        l5z = 0.1136;
        lfx = 0;
        lfz = 0.05716;

        l1m = 2.351;
        l2m = 1.318;
        l3m = 0.821;
        l4m = 0.769;
        l5m = 0.687;
        mg = 0.219;

        c1x = 0.01516;
        c1y = 0.00359;
        c1z = -0.03105;
        c2x = -0.01903;
        c2y = 0.0150;
        c2z = 0.11397;
        c3x = 0.00013;
        c3y = 0.02022;
        c3z = 0.10441;
        c4x = 0.00015;
        c4y = -0.02464;
        c4z = 0.05353;
        c5x = 0;
        c5y = 0.0012;
        c5z = 0.01648;
        cgx = 0;
        cgy = 0;
        cgz = 0.0289;

        I1x = 0.0029525;
        I1y = 0.0060091;
        I1z = 0.0058821;
        I2x = 0.0031145;
        I2y = 0.0005843;
        I2z = 0.0031631;
        I3x = 0.00172767;
        I3y = 0.00041967;
        I3z = 0.0018468;
        I4x = 0.0006764;
        I4y = 0.0010573;
        I4z = 0.0006610;
        I5x = 0.0001934;
        I5y = 0.0001602;
        I5z = 0.0000689;
        Igx = 0.0002324;
        Igy = 0.0003629;
        Igz = 0.0002067;
        sin00 = sin(pos[4])
        cos00 = cos(pos[4])
        cos01 = cos(pos[3])
        cos02 = cos(pos[1] - pos[2] + pos[3] - pos[4])
        sin01 = sin(pos[1] - pos[2] + pos[3] + pos[4])
        cos03 = cos(pos[1] - pos[2] + pos[3] + pos[4])
        sin02 = sin(pos[1] - pos[2] + pos[3] - pos[4])
        cos04 = cos(pos[1] - pos[2] + pos[3] - (2.0 * pos[4]))
        cos05 = cos(pos[1] - pos[2] + pos[3] + (2.0 * pos[4]))
        cos06 = cos(pos[2] - pos[3])
        sin03 = sin(pos[3])
        cos07 = cos(pos[1] - pos[2] + pos[3])
        sin04 = sin(pos[2] - pos[3])
        sin05 = sin(pos[1] - pos[2] - pos[4])
        cos08 = cos(pos[1] - pos[2] + pos[4])
        cos09 = cos(pos[1] - pos[2] - pos[4])
        sin06 = sin(pos[1] - pos[2] + pos[4])
        sin07 = sin(pos[1] - pos[2] + pos[3] - (2.0 * pos[4]))
        sin08 = sin(pos[1] - pos[2] + pos[3] + (2.0 * pos[4]))
        cos10 = cos(pos[1] - pos[2])
        cos11 = cos(pos[1] + pos[4])
        sin09 = sin(pos[1] - pos[4])
        sin10 = sin(pos[1] + pos[4])
        cos12 = cos(pos[1] - pos[4])
        sin11 = sin(pos[1] - pos[2] + pos[3])
        cos13 = cos(pos[1])
        sin12 = sin(pos[1] - pos[2])
        sin13 = sin(pos[1])
        sin14 = sin(pos[2])
        cos14 = cos(pos[2])
        M[0,0] = I1z + (I5y * pow(cos02 - cos03, 2.0) / 4.0) + (Igy * pow(cos02 - cos03, 2.0) / 4.0) + (
                    I3z * pow(cos10, 2.0)) + (I3x * pow(sin12, 2.0)) + (
                              l2m * pow(l2x + (c2x * cos13) - (c2z * sin13), 2.0)) + (
                              l5m * pow((c5x * cos02) + (c5y * sin02) - (c5x * cos03) + (c5y * sin01), 2.0) / 4.0) + (
                              mg * pow((cgx * cos02) + (cgy * sin02) - (cgx * cos03) + (cgy * sin01), 2.0) / 4.0) + (
                              l3m * pow(l2x - (l3z * sin13) + (c3x * cos10) - (c3z * sin12), 2.0)) + (l5m * pow(
            (c5z * cos02) + (l5z * cos02) + (l5x * sin02) + (2.0 * c5y * cos07) - (l4z * cos08) - (l3z * cos11) + (
                        l4z * cos09) - (2.0 * l2x * sin00) - (c5z * cos03) - (l5z * cos03) - (l5x * sin01) + (l3z * cos12),
            2.0) / 4.0) + (mg * pow(
            (cgz * cos02) + (l5z * cos02) + (l5x * sin02) + (2.0 * cgy * cos07) - (l4z * cos08) - (l3z * cos11) + (
                        l4z * cos09) - (2.0 * l2x * sin00) - (cgz * cos03) - (l5z * cos03) - (l5x * sin01) + (l3z * cos12),
            2.0) / 4.0) + (l4m * pow((c4z * sin11) - (c4x * cos07) - l2x + (l3z * sin13) + (l4z * sin12), 2.0)) + (
                              I4z * pow(cos07, 2.0)) + (I5z * pow(cos07, 2.0)) + (Igz * pow(cos07, 2.0)) + (
                              pow(c1x, 2.0) * l1m) + (pow(c1y, 2.0) * l1m) + (I4x * pow(sin11, 2.0)) + (l5m * pow(
            (c5z * sin02) - (l5x * cos02) + (l5z * sin02) - (2.0 * c5x * cos07) + (l4z * sin06) + (l3z * sin10) - (
                        2.0 * l2x * cos00) + (l4z * sin05) - (l5x * cos03) + (c5z * sin01) + (l5z * sin01) + (l3z * sin09),
            2.0) / 4.0) + (mg * pow(
            (cgz * sin02) - (l5x * cos02) + (l5z * sin02) - (2.0 * cgx * cos07) + (l4z * sin06) + (l3z * sin10) - (
                        2.0 * l2x * cos00) + (l4z * sin05) - (l5x * cos03) + (cgz * sin01) + (l5z * sin01) + (l3z * sin09),
            2.0) / 4.0) + (I5x * pow(sin01 + sin02, 2.0) / 4.0) + (Igx * pow(sin01 + sin02, 2.0) / 4.0) + (
                              I2z * pow(cos13, 2.0)) + (I2x * pow(sin13, 2.0)) + (
                              pow(c3y, 2.0) * l3m * pow(sin12, 2.0)) + (pow(c4y, 2.0) * l4m * pow(cos07, 2.0)) + (
                              pow(c4y, 2.0) * l4m * pow(sin11, 2.0)) + (pow(c2y, 2.0) * l2m * pow(cos13, 2.0)) + (
                              pow(c2y, 2.0) * l2m * pow(sin13, 2.0)) + (pow(c3y, 2.0) * l3m * pow(cos10, 2.0));
        M[0,1] = (I5x * cos05 / 4.0) - (I5x * cos04 / 4.0) + (I5y * cos04 / 4.0) - (I5y * cos05 / 4.0) - (
                    Igx * cos04 / 4.0) + (Igx * cos05 / 4.0) + (Igy * cos04 / 4.0) - (Igy * cos05 / 4.0) + (
                              pow(c5x, 2.0) * l5m * cos04 / 4.0) - (pow(c5x, 2.0) * l5m * cos05 / 4.0) - (
                              pow(c5y, 2.0) * l5m * cos04 / 4.0) + (pow(c5y, 2.0) * l5m * cos05 / 4.0) + (
                              pow(cgx, 2.0) * mg * cos04 / 4.0) - (pow(cgx, 2.0) * mg * cos05 / 4.0) - (
                              pow(cgy, 2.0) * mg * cos04 / 4.0) + (pow(cgy, 2.0) * mg * cos05 / 4.0) + (
                              c4y * c4z * l4m * cos07) + (c5y * l5m * l4z * cos08 / 2.0) + (
                              cgy * l4z * mg * cos08 / 2.0) + (c4x * c4y * l4m * sin11) + (
                              c5x * l5m * l4z * sin06 / 2.0) + (cgx * l4z * mg * sin06 / 2.0) + (
                              c5y * l5m * l3z * cos11 / 2.0) + (cgy * l3z * mg * cos11 / 2.0) + (
                              c5x * l5m * l3z * sin10 / 2.0) + (cgx * l3z * mg * sin10 / 2.0) + (
                              c2y * c2z * l2m * cos13) + (c3y * l3m * l3z * cos13) + (c4y * l4m * l3z * cos13) + (
                              c5y * l5m * l4z * cos09 / 2.0) + (c2x * c2y * l2m * sin13) + (
                              cgy * l4z * mg * cos09 / 2.0) + (c5y * c5z * l5m * cos03 / 2.0) + (
                              cgy * cgz * mg * cos03 / 2.0) - (c5x * l5m * l4z * sin05 / 2.0) - (
                              cgx * l4z * mg * sin05 / 2.0) - (c5x * l5m * l5x * cos03 / 2.0) + (
                              c5y * l5m * l5z * cos03 / 2.0) - (cgx * l5x * mg * cos03 / 2.0) + (
                              cgy * l5z * mg * cos03 / 2.0) + (c5x * c5z * l5m * sin01 / 2.0) + (
                              cgx * cgz * mg * sin01 / 2.0) + (c3y * c3z * l3m * cos10) + (
                              c5y * l5m * l5x * sin01 / 2.0) + (c5x * l5m * l5z * sin01 / 2.0) + (
                              cgy * l5x * mg * sin01 / 2.0) + (cgx * l5z * mg * sin01 / 2.0) + (c4y * l4m * l4z * cos10) + (
                              c5y * l5m * l3z * cos12 / 2.0) + (cgy * l3z * mg * cos12 / 2.0) + (
                              c3x * c3y * l3m * sin12) - (c5x * l5m * l3z * sin09 / 2.0) - (
                              cgx * l3z * mg * sin09 / 2.0) + (c5y * c5z * l5m * cos02 / 2.0) + (
                              cgy * cgz * mg * cos02 / 2.0) + (c5x * l5m * l5x * cos02 / 2.0) + (
                              c5y * l5m * l5z * cos02 / 2.0) + (cgx * l5x * mg * cos02 / 2.0) + (
                              cgy * l5z * mg * cos02 / 2.0) + (c5x * c5y * l5m * sin07 / 2.0) + (
                              c5x * c5y * l5m * sin08 / 2.0) - (c5x * c5z * l5m * sin02 / 2.0) + (
                              cgx * cgy * mg * sin07 / 2.0) + (cgx * cgy * mg * sin08 / 2.0) - (
                              cgx * cgz * mg * sin02 / 2.0) + (c5y * l5m * l5x * sin02 / 2.0) - (
                              c5x * l5m * l5z * sin02 / 2.0) + (cgy * l5x * mg * sin02 / 2.0) - (
                              cgx * l5z * mg * sin02 / 2.0)
        M[0,2] = (I5x * cos04 / 4.0) - (I5x * cos05 / 4.0) - (I5y * cos04 / 4.0) + (I5y * cos05 / 4.0) + (
                    Igx * cos04 / 4.0) - (Igx * cos05 / 4.0) - (Igy * cos04 / 4.0) + (Igy * cos05 / 4.0) - (
                              pow(c5x, 2.0) * l5m * cos04 / 4.0) + (pow(c5x, 2.0) * l5m * cos05 / 4.0) + (
                              pow(c5y, 2.0) * l5m * cos04 / 4.0) - (pow(c5y, 2.0) * l5m * cos05 / 4.0) - (
                              pow(cgx, 2.0) * mg * cos04 / 4.0) + (pow(cgx, 2.0) * mg * cos05 / 4.0) + (
                              pow(cgy, 2.0) * mg * cos04 / 4.0) - (pow(cgy, 2.0) * mg * cos05 / 4.0) - (
                              c4y * c4z * l4m * cos07) - (c5y * l5m * l4z * cos08 / 2.0) - (
                              cgy * l4z * mg * cos08 / 2.0) - (c4x * c4y * l4m * sin11) - (
                              c5x * l5m * l4z * sin06 / 2.0) - (cgx * l4z * mg * sin06 / 2.0) - (
                              c5y * l5m * l4z * cos09 / 2.0) - (cgy * l4z * mg * cos09 / 2.0) - (
                              c5y * c5z * l5m * cos03 / 2.0) - (cgy * cgz * mg * cos03 / 2.0) + (
                              c5x * l5m * l4z * sin05 / 2.0) + (cgx * l4z * mg * sin05 / 2.0) + (
                              c5x * l5m * l5x * cos03 / 2.0) - (c5y * l5m * l5z * cos03 / 2.0) + (
                              cgx * l5x * mg * cos03 / 2.0) - (cgy * l5z * mg * cos03 / 2.0) - (
                              c5x * c5z * l5m * sin01 / 2.0) - (cgx * cgz * mg * sin01 / 2.0) - (
                              c3y * c3z * l3m * cos10) - (c5y * l5m * l5x * sin01 / 2.0) - (
                              c5x * l5m * l5z * sin01 / 2.0) - (cgy * l5x * mg * sin01 / 2.0) - (
                              cgx * l5z * mg * sin01 / 2.0) - (c4y * l4m * l4z * cos10) - (c3x * c3y * l3m * sin12) - (
                              c5y * c5z * l5m * cos02 / 2.0) - (cgy * cgz * mg * cos02 / 2.0) - (
                              c5x * l5m * l5x * cos02 / 2.0) - (c5y * l5m * l5z * cos02 / 2.0) - (
                              cgx * l5x * mg * cos02 / 2.0) - (cgy * l5z * mg * cos02 / 2.0) - (
                              c5x * c5y * l5m * sin07 / 2.0) - (c5x * c5y * l5m * sin08 / 2.0) + (
                              c5x * c5z * l5m * sin02 / 2.0) - (cgx * cgy * mg * sin07 / 2.0) - (
                              cgx * cgy * mg * sin08 / 2.0) + (cgx * cgz * mg * sin02 / 2.0) - (
                              c5y * l5m * l5x * sin02 / 2.0) + (c5x * l5m * l5z * sin02 / 2.0) - (
                              cgy * l5x * mg * sin02 / 2.0) + (cgx * l5z * mg * sin02 / 2.0)
        M[0,3] = (I5x * cos05 / 4.0) - (I5x * cos04 / 4.0) + (I5y * cos04 / 4.0) - (I5y * cos05 / 4.0) - (
                    Igx * cos04 / 4.0) + (Igx * cos05 / 4.0) + (Igy * cos04 / 4.0) - (Igy * cos05 / 4.0) + (
                              pow(c5x, 2.0) * l5m * cos04 / 4.0) - (pow(c5x, 2.0) * l5m * cos05 / 4.0) - (
                              pow(c5y, 2.0) * l5m * cos04 / 4.0) + (pow(c5y, 2.0) * l5m * cos05 / 4.0) + (
                              pow(cgx, 2.0) * mg * cos04 / 4.0) - (pow(cgx, 2.0) * mg * cos05 / 4.0) - (
                              pow(cgy, 2.0) * mg * cos04 / 4.0) + (pow(cgy, 2.0) * mg * cos05 / 4.0) + (
                              c4y * c4z * l4m * cos07) + (c4x * c4y * l4m * sin11) + (c5y * c5z * l5m * cos03 / 2.0) + (
                              cgy * cgz * mg * cos03 / 2.0) - (c5x * l5m * l5x * cos03 / 2.0) + (
                              c5y * l5m * l5z * cos03 / 2.0) - (cgx * l5x * mg * cos03 / 2.0) + (
                              cgy * l5z * mg * cos03 / 2.0) + (c5x * c5z * l5m * sin01 / 2.0) + (
                              cgx * cgz * mg * sin01 / 2.0) + (c5y * l5m * l5x * sin01 / 2.0) + (
                              c5x * l5m * l5z * sin01 / 2.0) + (cgy * l5x * mg * sin01 / 2.0) + (
                              cgx * l5z * mg * sin01 / 2.0) + (c5y * c5z * l5m * cos02 / 2.0) + (
                              cgy * cgz * mg * cos02 / 2.0) + (c5x * l5m * l5x * cos02 / 2.0) + (
                              c5y * l5m * l5z * cos02 / 2.0) + (cgx * l5x * mg * cos02 / 2.0) + (
                              cgy * l5z * mg * cos02 / 2.0) + (c5x * c5y * l5m * sin07 / 2.0) + (
                              c5x * c5y * l5m * sin08 / 2.0) - (c5x * c5z * l5m * sin02 / 2.0) + (
                              cgx * cgy * mg * sin07 / 2.0) + (cgx * cgy * mg * sin08 / 2.0) - (
                              cgx * cgz * mg * sin02 / 2.0) + (c5y * l5m * l5x * sin02 / 2.0) - (
                              c5x * l5m * l5z * sin02 / 2.0) + (cgy * l5x * mg * sin02 / 2.0) - (
                              cgx * l5z * mg * sin02 / 2.0);
        M[0,4] = (I5z * cos07) + (Igz * cos07) + (pow(c5x, 2.0) * l5m * cos07) + (pow(c5y, 2.0) * l5m * cos07) + (
                    pow(cgx, 2.0) * mg * cos07) + (pow(cgy, 2.0) * mg * cos07) - (c5y * l5m * l4z * cos08 / 2.0) - (
                              cgy * l4z * mg * cos08 / 2.0) - (c5x * l5m * l4z * sin06 / 2.0) - (
                              cgx * l4z * mg * sin06 / 2.0) - (c5y * l5m * l3z * cos11 / 2.0) - (
                              cgy * l3z * mg * cos11 / 2.0) - (c5x * l5m * l3z * sin10 / 2.0) - (
                              cgx * l3z * mg * sin10 / 2.0) + (c5x * l5m * l2x * cos00) + (cgx * l2x * mg * cos00) + (
                              c5y * l5m * l4z * cos09 / 2.0) + (cgy * l4z * mg * cos09 / 2.0) - (
                              c5y * l5m * l2x * sin00) - (cgy * l2x * mg * sin00) - (c5y * c5z * l5m * cos03 / 2.0) - (
                              cgy * cgz * mg * cos03 / 2.0) - (c5x * l5m * l4z * sin05 / 2.0) - (
                              cgx * l4z * mg * sin05 / 2.0) + (c5x * l5m * l5x * cos03 / 2.0) - (
                              c5y * l5m * l5z * cos03 / 2.0) + (cgx * l5x * mg * cos03 / 2.0) - (
                              cgy * l5z * mg * cos03 / 2.0) - (c5x * c5z * l5m * sin01 / 2.0) - (
                              cgx * cgz * mg * sin01 / 2.0) - (c5y * l5m * l5x * sin01 / 2.0) - (
                              c5x * l5m * l5z * sin01 / 2.0) - (cgy * l5x * mg * sin01 / 2.0) - (
                              cgx * l5z * mg * sin01 / 2.0) + (c5y * l5m * l3z * cos12 / 2.0) + (
                              cgy * l3z * mg * cos12 / 2.0) - (c5x * l5m * l3z * sin09 / 2.0) - (
                              cgx * l3z * mg * sin09 / 2.0) + (c5y * c5z * l5m * cos02 / 2.0) + (
                              cgy * cgz * mg * cos02 / 2.0) + (c5x * l5m * l5x * cos02 / 2.0) + (
                              c5y * l5m * l5z * cos02 / 2.0) + (cgx * l5x * mg * cos02 / 2.0) + (
                              cgy * l5z * mg * cos02 / 2.0) - (c5x * c5z * l5m * sin02 / 2.0) - (
                              cgx * cgz * mg * sin02 / 2.0) + (c5y * l5m * l5x * sin02 / 2.0) - (
                              c5x * l5m * l5z * sin02 / 2.0) + (cgy * l5x * mg * sin02 / 2.0) - (
                              cgx * l5z * mg * sin02 / 2.0);
        M[1,0] = M[0,1]
        M[1,1] = I2y + I3y + I4y + (l4m * pow(c4z + (l4z * cos01) + (l3z * cos06), 2.0)) + (
                    l4m * pow(c4x + (l4z * sin03) - (l3z * sin04), 2.0)) + (pow(c2x, 2.0) * l2m) + (
                              pow(c2z, 2.0) * l2m) + (l3m * pow(c3z + (l3z * cos14), 2.0)) + (
                              l3m * pow(c3x - (l3z * sin14), 2.0)) + (
                              l5m * pow(l5x + (c5x * cos00) - (c5y * sin00) + (l4z * sin03) - (l3z * sin04), 2.0)) + (
                              mg * pow(l5x + (cgx * cos00) - (cgy * sin00) + (l4z * sin03) - (l3z * sin04), 2.0)) + (
                              I5y * pow(cos00, 2.0)) + (Igy * pow(cos00, 2.0)) + (I5x * pow(sin00, 2.0)) + (
                              Igx * pow(sin00, 2.0)) + (
                              l5m * pow(cos00, 2.0) * pow(c5z + l5z + (l4z * cos01) + (l3z * cos06), 2.0)) + (
                              mg * pow(cos00, 2.0) * pow(cgz + l5z + (l4z * cos01) + (l3z * cos06), 2.0)) + (
                              l5m * pow(sin00, 2.0) * pow(c5z + l5z + (l4z * cos01) + (l3z * cos06), 2.0)) + (
                              mg * pow(sin00, 2.0) * pow(cgz + l5z + (l4z * cos01) + (l3z * cos06), 2.0))
        M[1,2] = -I3y - I4y - (I5y * pow(cos00, 2.0)) - (Igy * pow(cos00, 2.0)) - (I5x * pow(sin00, 2.0)) - (
                    Igx * pow(sin00, 2.0)) - (l5m * (l5x + (c5x * cos00) - (c5y * sin00) + (l4z * sin03)) * (
                    l5x + (c5x * cos00) - (c5y * sin00) + (l4z * sin03) - (l3z * sin04))) - (
                              mg * (l5x + (cgx * cos00) - (cgy * sin00) + (l4z * sin03)) * (
                                  l5x + (cgx * cos00) - (cgy * sin00) + (l4z * sin03) - (l3z * sin04))) - (
                              l4m * (c4z + (l4z * cos01)) * (c4z + (l4z * cos01) + (l3z * cos06))) - (
                              c3z * l3m * (c3z + (l3z * cos14))) - (
                              l4m * (c4x + (l4z * sin03)) * (c4x + (l4z * sin03) - (l3z * sin04))) - (
                              c3x * l3m * (c3x - (l3z * sin14))) - (
                              l5m * pow(cos00, 2.0) * (c5z + l5z + (l4z * cos01)) * (
                                  c5z + l5z + (l4z * cos01) + (l3z * cos06))) - (
                              mg * pow(cos00, 2.0) * (cgz + l5z + (l4z * cos01)) * (
                                  cgz + l5z + (l4z * cos01) + (l3z * cos06))) - (
                              l5m * pow(sin00, 2.0) * (c5z + l5z + (l4z * cos01)) * (
                                  c5z + l5z + (l4z * cos01) + (l3z * cos06))) - (
                              mg * pow(sin00, 2.0) * (cgz + l5z + (l4z * cos01)) * (
                                  cgz + l5z + (l4z * cos01) + (l3z * cos06)))
        M[1,3] = I4y + (I5y * pow(cos00, 2.0)) + (Igy * pow(cos00, 2.0)) + (I5x * pow(sin00, 2.0)) + (
                    Igx * pow(sin00, 2.0)) + (c4z * l4m * (c4z + (l4z * cos01) + (l3z * cos06))) + (
                              l5m * (l5x + (c5x * cos00) - (c5y * sin00)) * (
                                  l5x + (c5x * cos00) - (c5y * sin00) + (l4z * sin03) - (l3z * sin04))) + (
                              mg * (l5x + (cgx * cos00) - (cgy * sin00)) * (
                                  l5x + (cgx * cos00) - (cgy * sin00) + (l4z * sin03) - (l3z * sin04))) + (
                              c4x * l4m * (c4x + (l4z * sin03) - (l3z * sin04))) + (
                              l5m * pow(cos00, 2.0) * (c5z + l5z) * (c5z + l5z + (l4z * cos01) + (l3z * cos06))) + (
                              mg * pow(cos00, 2.0) * (cgz + l5z) * (cgz + l5z + (l4z * cos01) + (l3z * cos06))) + (
                              l5m * pow(sin00, 2.0) * (c5z + l5z) * (c5z + l5z + (l4z * cos01) + (l3z * cos06))) + (
                              mg * pow(sin00, 2.0) * (cgz + l5z) * (cgz + l5z + (l4z * cos01) + (l3z * cos06)))
        M[1,4] = (c5y * l5m * cos00 * (c5z + l5z + (l4z * cos01) + (l3z * cos06))) + (
                    cgy * mg * cos00 * (cgz + l5z + (l4z * cos01) + (l3z * cos06))) + (
                              c5x * l5m * sin00 * (c5z + l5z + (l4z * cos01) + (l3z * cos06))) + (
                              cgx * mg * sin00 * (cgz + l5z + (l4z * cos01) + (l3z * cos06)))
        M[2,0] = M[0,2]
        M[2,1] = M[1,2]
        M[2,2] = I3y + I4y + (pow(c3x, 2.0) * l3m) + (pow(c3z, 2.0) * l3m) + (
                    l4m * pow(c4z + (l4z * cos01), 2.0)) + (l4m * pow(c4x + (l4z * sin03), 2.0)) + (
                              I5y * pow(cos00, 2.0)) + (Igy * pow(cos00, 2.0)) + (
                              l5m * pow(l5x + (c5x * cos00) - (c5y * sin00) + (l4z * sin03), 2.0)) + (
                              mg * pow(l5x + (cgx * cos00) - (cgy * sin00) + (l4z * sin03), 2.0)) + (
                              I5x * pow(sin00, 2.0)) + (Igx * pow(sin00, 2.0)) + (
                              l5m * pow(cos00, 2.0) * pow(c5z + l5z + (l4z * cos01), 2.0)) + (
                              mg * pow(cos00, 2.0) * pow(cgz + l5z + (l4z * cos01), 2.0)) + (
                              l5m * pow(sin00, 2.0) * pow(c5z + l5z + (l4z * cos01), 2.0)) + (
                              mg * pow(sin00, 2.0) * pow(cgz + l5z + (l4z * cos01), 2.0));
        M[2,3] = -I4y - (I5y * pow(cos00, 2.0)) - (Igy * pow(cos00, 2.0)) - (I5x * pow(sin00, 2.0)) - (
                    Igx * pow(sin00, 2.0)) - (c4z * l4m * (c4z + (l4z * cos01))) - (
                              l5m * (l5x + (c5x * cos00) - (c5y * sin00)) * (
                                  l5x + (c5x * cos00) - (c5y * sin00) + (l4z * sin03))) - (
                              mg * (l5x + (cgx * cos00) - (cgy * sin00)) * (
                                  l5x + (cgx * cos00) - (cgy * sin00) + (l4z * sin03))) - (
                              c4x * l4m * (c4x + (l4z * sin03))) - (
                              l5m * pow(cos00, 2.0) * (c5z + l5z) * (c5z + l5z + (l4z * cos01))) - (
                              mg * pow(cos00, 2.0) * (cgz + l5z) * (cgz + l5z + (l4z * cos01))) - (
                              l5m * pow(sin00, 2.0) * (c5z + l5z) * (c5z + l5z + (l4z * cos01))) - (
                              mg * pow(sin00, 2.0) * (cgz + l5z) * (cgz + l5z + (l4z * cos01)))
        M[2,4] = -(c5y * l5m * cos00 * (c5z + l5z + (l4z * cos01))) - (cgy * mg * cos00 * (cgz + l5z + (l4z * cos01))) - (
                    c5x * l5m * sin00 * (c5z + l5z + (l4z * cos01))) - (cgx * mg * sin00 * (cgz + l5z + (l4z * cos01)))
        M[3,0] = M[0,3]
        M[3,1] = M[1,3]
        M[3,2] = M[2,3]
        M[3,3] = I4y + (l5m * pow(l5x + (c5x * cos00) - (c5y * sin00), 2.0)) + (
                    mg * pow(l5x + (cgx * cos00) - (cgy * sin00), 2.0)) + (pow(c4x, 2.0) * l4m) + (
                              pow(c4z, 2.0) * l4m) + (I5y * pow(cos00, 2.0)) + (Igy * pow(cos00, 2.0)) + (
                              I5x * pow(sin00, 2.0)) + (Igx * pow(sin00, 2.0)) + (
                              l5m * pow(cos00, 2.0) * pow(c5z + l5z, 2.0)) + (
                              mg * pow(cos00, 2.0) * pow(cgz + l5z, 2.0)) + (
                              l5m * pow(sin00, 2.0) * pow(c5z + l5z, 2.0)) + (
                              mg * pow(sin00, 2.0) * pow(cgz + l5z, 2.0))
        M[3,4] = (c5y * l5m * cos00 * (c5z + l5z)) + (cgy * mg * cos00 * (cgz + l5z)) + (
                    c5x * l5m * sin00 * (c5z + l5z)) + (cgx * mg * sin00 * (cgz + l5z))
        M[4,0] = M[0,4]
        M[4,1] = M[1,4]
        M[4,2] = M[2,4]
        M[4,3] = M[3,4]
        M[4,4]  = I5z + Igz + (pow(c5x, 2.0) * l5m) + (pow(c5y, 2.0) * l5m) + (pow(cgx, 2.0) * mg) + (
                    pow(cgy, 2.0) * mg)


    def C_matrix(self):
        pos = self.pos
        vel = self.vel
        C = self.C
        lbx = 0;
        lbz = 0;
        l1x = 0.024;
        l1z = 0.115;
        l2x = 0.033;
        l2z = 0;
        l3x = 0;
        l3z = 0.155;
        l4x = 0;
        l4z = 0.135;
        l5x = 0;
        l5z = 0.1136;
        lfx = 0;
        lfz = 0.05716;

        l1m = 2.351;
        l2m = 1.318;
        l3m = 0.821;
        l4m = 0.769;
        l5m = 0.687;
        mg = 0.219;

        c1x = 0.01516;
        c1y = 0.00359;
        c1z = -0.03105;
        c2x = -0.01903;
        c2y = 0.0150;
        c2z = 0.11397;
        c3x = 0.00013;
        c3y = 0.02022;
        c3z = 0.10441;
        c4x = 0.00015;
        c4y = -0.02464;
        c4z = 0.05353;
        c5x = 0;
        c5y = 0.0012;
        c5z = 0.01648;
        cgx = 0;
        cgy = 0;
        cgz = 0.0289;

        I1x = 0.0029525;
        I1y = 0.0060091;
        I1z = 0.0058821;
        I2x = 0.0031145;
        I2y = 0.0005843;
        I2z = 0.0031631;
        I3x = 0.00172767;
        I3y = 0.00041967;
        I3z = 0.0018468;
        I4x = 0.0006764;
        I4y = 0.0010573;
        I4z = 0.0006610;
        I5x = 0.0001934;
        I5y = 0.0001602;
        I5z = 0.0000689;
        Igx = 0.0002324;
        Igy = 0.0003629;
        Igz = 0.0002067;
        sin00 = sin(2.0 * pos[4])
        sin01 = sin(pos[1] - pos[2] + pos[3] - (2.0 * pos[4]))
        sin02 = sin(pos[1] - pos[2] + pos[3] + (2.0 * pos[4]))
        cos00 = cos(pos[1] - pos[2] + pos[3] - pos[4])
        sin03 = sin(pos[1] - pos[2] + pos[3] + pos[4])
        cos01 = cos(pos[1] - pos[2] + pos[3] + pos[4])
        sin04 = sin(pos[1] - pos[2] + pos[3] - pos[4])
        sin05 = sin((2.0 * pos[1]) - (2.0 * pos[2]) + (2.0 * pos[3]))
        sin06 = sin(pos[1] - pos[2] + pos[3])
        cos02 = cos(pos[4])
        sin07 = sin(pos[4])
        sin08 = sin(pos[3])
        sin09 = sin((2.0 * pos[1]) - (2.0 * pos[2]) + (2.0 * pos[3]) + (2.0 * pos[4]))
        sin10 = sin((2.0 * pos[1]) - (2.0 * pos[2]) + (2.0 * pos[3]) - (2.0 * pos[4]))
        cos03 = cos(pos[3])
        sin11 = sin(pos[2] - pos[3])
        sin12 = sin((2.0 * pos[1]) - (2.0 * pos[2]) + (2.0 * pos[3]) - pos[4])
        sin13 = sin((2.0 * pos[1]) - (2.0 * pos[2]) + (2.0 * pos[3]) + pos[4])
        cos04 = cos(pos[1] - pos[2] + pos[3])
        cos05 = cos((2.0 * pos[1]) - (2.0 * pos[2]) + (2.0 * pos[3]) - pos[4])
        cos06 = cos((2.0 * pos[1]) - (2.0 * pos[2]) + (2.0 * pos[3]) + pos[4])
        cos07 = cos(pos[3] + pos[4])
        cos08 = cos(pos[1] - pos[2] + pos[3] - (2.0 * pos[4]))
        cos09 = cos(pos[3] - pos[4])
        sin14 = sin(pos[3] + pos[4])
        sin15 = sin(pos[3] - pos[4])
        cos10 = cos(pos[1] - pos[2] + pos[3] + (2.0 * pos[4]))
        cos11 = cos(2.0 * pos[4])
        cos12 = cos(pos[2] - pos[3])
        sin16 = sin(pos[2] - pos[3] + pos[4])
        cos13 = cos(pos[2] - pos[3] - pos[4])
        cos14 = cos(pos[2] - pos[3] + pos[4])
        sin17 = sin(pos[2] - pos[3] - pos[4])
        sin18 = sin((2.0 * pos[1]) - (2.0 * pos[2]))
        sin19 = sin((2.0 * pos[1]) - pos[2] + pos[3])
        cos15 = cos((2.0 * pos[1]) - (2.0 * pos[2]) + (2.0 * pos[3]))
        sin20 = sin((2.0 * pos[1]) - (2.0 * pos[2]) + pos[3])
        cos16 = cos(pos[1] - pos[2] + pos[4])
        sin21 = sin(pos[1] - pos[2] + pos[4])
        sin22 = sin(pos[1] - pos[2] - pos[4])
        cos17 = cos(pos[1] - pos[2] - pos[4])
        cos18 = cos(pos[1] - pos[2])
        sin23 = sin(2.0 * pos[1])
        sin24 = sin(pos[2])
        cos19 = cos((2.0 * pos[1]) - (2.0 * pos[2]) + pos[3])
        cos20 = cos((2.0 * pos[1]) - pos[2] + pos[3])
        cos21 = cos((2.0 * pos[1]) - (2.0 * pos[2]) + pos[3] - pos[4])
        sin25 = sin((2.0 * pos[1]) - (2.0 * pos[2]) + pos[3] - pos[4])
        sin26 = sin((2.0 * pos[1]) - pos[2] + pos[3] + pos[4])
        cos22 = cos((2.0 * pos[1]) - pos[2] + pos[3] - pos[4])
        cos23 = cos((2.0 * pos[1]) - (2.0 * pos[2]) + (2.0 * pos[3]) - (2.0 * pos[4]))
        sin27 = sin((2.0 * pos[1]) - pos[2] + pos[3] - pos[4])
        sin28 = sin((2.0 * pos[1]) - (2.0 * pos[2]) + pos[3] + pos[4])
        cos24 = cos((2.0 * pos[1]) - (2.0 * pos[2]) + pos[3] + pos[4])
        cos25 = cos((2.0 * pos[1]) - pos[2] + pos[3] + pos[4])
        cos26 = cos((2.0 * pos[1]) - (2.0 * pos[2]) + (2.0 * pos[3]) + (2.0 * pos[4]))
        sin29 = sin((2.0 * pos[1]) - pos[2])
        sin30 = sin(pos[1] - pos[2])
        sin31 = sin(pos[1] + pos[4])
        cos27 = cos(pos[1])
        cos28 = cos(pos[1] + pos[4])
        sin32 = sin(pos[1] - pos[4])
        cos29 = cos(pos[1] - pos[4])
        sin33 = sin(pos[1])
        cos30 = cos(pos[2])
        cos31 = cos((2.0 * pos[1]) - (2.0 * pos[2]))
        cos32 = cos((2.0 * pos[1]) - pos[2])
        cos33 = cos(2.0 * pos[1])

        C[0, 0] = (I5x * vel[1] * sin10 / 8.0) + (I5x * vel[1] * sin09 / 8.0) - (I5x * vel[2] * sin10 / 8.0) - (
                    I5x * vel[2] * sin09 / 8.0) + (I5x * vel[3] * sin10 / 8.0) + (I5x * vel[3] * sin09 / 8.0) - (
                              I5x * vel[4] * sin10 / 8.0) + (I5x * vel[4] * sin09 / 8.0) - (
                              I5y * vel[1] * sin10 / 8.0) - (I5y * vel[1] * sin09 / 8.0) + (
                              I5y * vel[2] * sin10 / 8.0) + (I5y * vel[2] * sin09 / 8.0) - (
                              I5y * vel[3] * sin10 / 8.0) - (I5y * vel[3] * sin09 / 8.0) + (
                              I5y * vel[4] * sin10 / 8.0) - (I5y * vel[4] * sin09 / 8.0) + (
                              Igx * vel[1] * sin10 / 8.0) + (Igx * vel[1] * sin09 / 8.0) - (
                              Igx * vel[2] * sin10 / 8.0) - (Igx * vel[2] * sin09 / 8.0) + (
                              Igx * vel[3] * sin10 / 8.0) + (Igx * vel[3] * sin09 / 8.0) - (
                              Igx * vel[4] * sin10 / 8.0) + (Igx * vel[4] * sin09 / 8.0) - (
                              Igy * vel[1] * sin10 / 8.0) - (Igy * vel[1] * sin09 / 8.0) + (
                              Igy * vel[2] * sin10 / 8.0) + (Igy * vel[2] * sin09 / 8.0) - (
                              Igy * vel[3] * sin10 / 8.0) - (Igy * vel[3] * sin09 / 8.0) + (
                              Igy * vel[4] * sin10 / 8.0) - (Igy * vel[4] * sin09 / 8.0) + (
                              I2x * vel[1] * sin23 / 2.0) - (I5x * vel[4] * sin00 / 4.0) + (
                              I5y * vel[4] * sin00 / 4.0) - (I2z * vel[1] * sin23 / 2.0) - (
                              Igx * vel[4] * sin00 / 4.0) + (Igy * vel[4] * sin00 / 4.0) + (
                              I4x * vel[1] * sin05 / 2.0) - (I4x * vel[2] * sin05 / 2.0) + (
                              I5x * vel[1] * sin05 / 4.0) + (I4x * vel[3] * sin05 / 2.0) - (
                              I5x * vel[2] * sin05 / 4.0) + (I5x * vel[3] * sin05 / 4.0) + (
                              I5y * vel[1] * sin05 / 4.0) - (I5y * vel[2] * sin05 / 4.0) + (
                              I5y * vel[3] * sin05 / 4.0) - (I4z * vel[1] * sin05 / 2.0) + (
                              I4z * vel[2] * sin05 / 2.0) - (I5z * vel[1] * sin05 / 2.0) - (
                              I4z * vel[3] * sin05 / 2.0) + (I5z * vel[2] * sin05 / 2.0) - (
                              I5z * vel[3] * sin05 / 2.0) + (Igx * vel[1] * sin05 / 4.0) - (
                              Igx * vel[2] * sin05 / 4.0) + (Igx * vel[3] * sin05 / 4.0) + (
                              Igy * vel[1] * sin05 / 4.0) - (Igy * vel[2] * sin05 / 4.0) + (
                              Igy * vel[3] * sin05 / 4.0) - (Igz * vel[1] * sin05 / 2.0) + (
                              Igz * vel[2] * sin05 / 2.0) - (Igz * vel[3] * sin05 / 2.0) + (
                              I3x * vel[1] * sin18 / 2.0) - (I3x * vel[2] * sin18 / 2.0) - (
                              I3z * vel[1] * sin18 / 2.0) + (I3z * vel[2] * sin18 / 2.0) - (
                              pow(c5x, 2.0) * l5m * vel[1] * sin10 / 8.0) - (
                              pow(c5x, 2.0) * l5m * vel[1] * sin09 / 8.0) + (
                              pow(c5x, 2.0) * l5m * vel[2] * sin10 / 8.0) + (
                              pow(c5x, 2.0) * l5m * vel[2] * sin09 / 8.0) - (
                              pow(c5x, 2.0) * l5m * vel[3] * sin10 / 8.0) - (
                              pow(c5x, 2.0) * l5m * vel[3] * sin09 / 8.0) + (
                              pow(c5x, 2.0) * l5m * vel[4] * sin10 / 8.0) - (
                              pow(c5x, 2.0) * l5m * vel[4] * sin09 / 8.0) + (
                              pow(c5y, 2.0) * l5m * vel[1] * sin10 / 8.0) + (
                              pow(c5y, 2.0) * l5m * vel[1] * sin09 / 8.0) - (
                              pow(c5y, 2.0) * l5m * vel[2] * sin10 / 8.0) - (
                              pow(c5y, 2.0) * l5m * vel[2] * sin09 / 8.0) + (
                              pow(c5y, 2.0) * l5m * vel[3] * sin10 / 8.0) + (
                              pow(c5y, 2.0) * l5m * vel[3] * sin09 / 8.0) - (
                              pow(c5y, 2.0) * l5m * vel[4] * sin10 / 8.0) + (
                              pow(c5y, 2.0) * l5m * vel[4] * sin09 / 8.0) - (
                              pow(cgx, 2.0) * mg * vel[1] * sin10 / 8.0) - (
                              pow(cgx, 2.0) * mg * vel[1] * sin09 / 8.0) + (
                              pow(cgx, 2.0) * mg * vel[2] * sin10 / 8.0) + (
                              pow(cgx, 2.0) * mg * vel[2] * sin09 / 8.0) - (
                              pow(cgx, 2.0) * mg * vel[3] * sin10 / 8.0) - (
                              pow(cgx, 2.0) * mg * vel[3] * sin09 / 8.0) + (
                              pow(cgx, 2.0) * mg * vel[4] * sin10 / 8.0) - (
                              pow(cgx, 2.0) * mg * vel[4] * sin09 / 8.0) + (
                              pow(cgy, 2.0) * mg * vel[1] * sin10 / 8.0) + (
                              pow(cgy, 2.0) * mg * vel[1] * sin09 / 8.0) - (
                              pow(cgy, 2.0) * mg * vel[2] * sin10 / 8.0) - (
                              pow(cgy, 2.0) * mg * vel[2] * sin09 / 8.0) + (
                              pow(cgy, 2.0) * mg * vel[3] * sin10 / 8.0) + (
                              pow(cgy, 2.0) * mg * vel[3] * sin09 / 8.0) - (
                              pow(cgy, 2.0) * mg * vel[4] * sin10 / 8.0) + (
                              pow(cgy, 2.0) * mg * vel[4] * sin09 / 8.0) - (
                              pow(c2x, 2.0) * l2m * vel[1] * sin23 / 2.0) + (
                              pow(c5x, 2.0) * l5m * vel[4] * sin00 / 4.0) - (
                              pow(c5y, 2.0) * l5m * vel[4] * sin00 / 4.0) + (
                              pow(c2z, 2.0) * l2m * vel[1] * sin23 / 2.0) + (
                              pow(cgx, 2.0) * mg * vel[4] * sin00 / 4.0) - (
                              pow(cgy, 2.0) * mg * vel[4] * sin00 / 4.0) - (
                              pow(c4x, 2.0) * l4m * vel[1] * sin05 / 2.0) + (
                              pow(c4x, 2.0) * l4m * vel[2] * sin05 / 2.0) - (
                              pow(c4x, 2.0) * l4m * vel[3] * sin05 / 2.0) - (
                              pow(c5x, 2.0) * l5m * vel[1] * sin05 / 4.0) + (
                              pow(c5x, 2.0) * l5m * vel[2] * sin05 / 4.0) - (
                              pow(c5x, 2.0) * l5m * vel[3] * sin05 / 4.0) - (
                              pow(c5y, 2.0) * l5m * vel[1] * sin05 / 4.0) + (
                              pow(c5y, 2.0) * l5m * vel[2] * sin05 / 4.0) - (
                              pow(c5y, 2.0) * l5m * vel[3] * sin05 / 4.0) + (
                              pow(c4z, 2.0) * l4m * vel[1] * sin05 / 2.0) - (
                              pow(c4z, 2.0) * l4m * vel[2] * sin05 / 2.0) + (
                              pow(c4z, 2.0) * l4m * vel[3] * sin05 / 2.0) + (
                              pow(c5z, 2.0) * l5m * vel[1] * sin05 / 2.0) - (
                              pow(c5z, 2.0) * l5m * vel[2] * sin05 / 2.0) + (
                              pow(c5z, 2.0) * l5m * vel[3] * sin05 / 2.0) - (
                              pow(cgx, 2.0) * mg * vel[1] * sin05 / 4.0) + (
                              pow(cgx, 2.0) * mg * vel[2] * sin05 / 4.0) - (
                              pow(cgx, 2.0) * mg * vel[3] * sin05 / 4.0) - (
                              pow(cgy, 2.0) * mg * vel[1] * sin05 / 4.0) + (
                              pow(cgy, 2.0) * mg * vel[2] * sin05 / 4.0) - (
                              pow(cgy, 2.0) * mg * vel[3] * sin05 / 4.0) + (
                              pow(cgz, 2.0) * mg * vel[1] * sin05 / 2.0) - (
                              pow(cgz, 2.0) * mg * vel[2] * sin05 / 2.0) + (
                              pow(cgz, 2.0) * mg * vel[3] * sin05 / 2.0) + (
                              l3m * pow(l3z, 2.0) * vel[1] * sin23 / 2.0) + (
                              l4m * pow(l3z, 2.0) * vel[1] * sin23 / 2.0) + (
                              l5m * pow(l3z, 2.0) * vel[1] * sin23 / 2.0) + (
                              pow(l3z, 2.0) * mg * vel[1] * sin23 / 2.0) - (
                              l5m * pow(l5x, 2.0) * vel[1] * sin05 / 2.0) + (
                              l5m * pow(l5x, 2.0) * vel[2] * sin05 / 2.0) - (
                              l5m * pow(l5x, 2.0) * vel[3] * sin05 / 2.0) + (
                              l5m * pow(l5z, 2.0) * vel[1] * sin05 / 2.0) - (
                              l5m * pow(l5z, 2.0) * vel[2] * sin05 / 2.0) + (
                              l5m * pow(l5z, 2.0) * vel[3] * sin05 / 2.0) - (
                              pow(l5x, 2.0) * mg * vel[1] * sin05 / 2.0) + (
                              pow(l5x, 2.0) * mg * vel[2] * sin05 / 2.0) - (
                              pow(l5x, 2.0) * mg * vel[3] * sin05 / 2.0) + (
                              pow(l5z, 2.0) * mg * vel[1] * sin05 / 2.0) - (
                              pow(l5z, 2.0) * mg * vel[2] * sin05 / 2.0) + (
                              pow(l5z, 2.0) * mg * vel[3] * sin05 / 2.0) - (
                              pow(c3x, 2.0) * l3m * vel[1] * sin18 / 2.0) + (
                              pow(c3x, 2.0) * l3m * vel[2] * sin18 / 2.0) + (
                              pow(c3z, 2.0) * l3m * vel[1] * sin18 / 2.0) - (
                              pow(c3z, 2.0) * l3m * vel[2] * sin18 / 2.0) + (
                              l4m * pow(l4z, 2.0) * vel[1] * sin18 / 2.0) - (
                              l4m * pow(l4z, 2.0) * vel[2] * sin18 / 2.0) + (
                              l5m * pow(l4z, 2.0) * vel[1] * sin18 / 2.0) - (
                              l5m * pow(l4z, 2.0) * vel[2] * sin18 / 2.0) + (
                              pow(l4z, 2.0) * mg * vel[1] * sin18 / 2.0) - (
                              pow(l4z, 2.0) * mg * vel[2] * sin18 / 2.0) - (l3z * l4z * mg * vel[2] * sin24 / 2.0) - (
                              l4z * l5z * mg * vel[3] * sin08 / 2.0) - (c5y * l5m * l2x * vel[1] * cos01 / 2.0) + (
                              c5y * l5m * l2x * vel[2] * cos01 / 2.0) - (c5y * l5m * l2x * vel[3] * cos01 / 2.0) - (
                              c5y * l5m * l2x * vel[4] * cos01 / 2.0) - (cgy * l2x * mg * vel[1] * cos01 / 2.0) + (
                              cgy * l2x * mg * vel[2] * cos01 / 2.0) - (cgy * l2x * mg * vel[3] * cos01 / 2.0) - (
                              cgy * l2x * mg * vel[4] * cos01 / 2.0) + (l5m * l3z * l5z * vel[1] * sin19) - (
                              l5m * l3z * l5z * vel[2] * sin19 / 2.0) + (l5m * l3z * l5z * vel[3] * sin19 / 2.0) + (
                              l5m * l4z * l5z * vel[1] * sin20) - (l5m * l4z * l5z * vel[2] * sin20) + (
                              l5m * l4z * l5z * vel[3] * sin20 / 2.0) + (l3z * l5z * mg * vel[1] * sin19) - (
                              l3z * l5z * mg * vel[2] * sin19 / 2.0) + (l3z * l5z * mg * vel[3] * sin19 / 2.0) + (
                              l4z * l5z * mg * vel[1] * sin20) - (l4z * l5z * mg * vel[2] * sin20) + (
                              l4z * l5z * mg * vel[3] * sin20 / 2.0) - (c5x * l5m * l2x * vel[1] * sin03 / 2.0) + (
                              c5x * l5m * l2x * vel[2] * sin03 / 2.0) - (c5x * l5m * l2x * vel[3] * sin03 / 2.0) - (
                              c5x * l5m * l2x * vel[4] * sin03 / 2.0) - (cgx * l2x * mg * vel[1] * sin03 / 2.0) + (
                              cgx * l2x * mg * vel[2] * sin03 / 2.0) - (cgx * l2x * mg * vel[3] * sin03 / 2.0) - (
                              cgx * l2x * mg * vel[4] * sin03 / 2.0) + (c5x * c5y * l5m * vel[1] * cos23 / 4.0) - (
                              c5x * c5y * l5m * vel[1] * cos26 / 4.0) - (c5x * c5y * l5m * vel[2] * cos23 / 4.0) + (
                              c5x * c5y * l5m * vel[2] * cos26 / 4.0) + (c5x * c5y * l5m * vel[3] * cos23 / 4.0) - (
                              c5x * c5y * l5m * vel[3] * cos26 / 4.0) - (c5x * c5y * l5m * vel[4] * cos23 / 4.0) - (
                              c5x * c5y * l5m * vel[4] * cos26 / 4.0) - (c5x * c5z * l5m * vel[1] * cos05 / 2.0) + (
                              c5x * c5z * l5m * vel[2] * cos05 / 2.0) - (c5x * c5z * l5m * vel[3] * cos05 / 2.0) + (
                              c5x * c5z * l5m * vel[4] * cos05 / 4.0) + (cgx * cgy * mg * vel[1] * cos23 / 4.0) - (
                              cgx * cgy * mg * vel[1] * cos26 / 4.0) - (cgx * cgy * mg * vel[2] * cos23 / 4.0) + (
                              cgx * cgy * mg * vel[2] * cos26 / 4.0) + (cgx * cgy * mg * vel[3] * cos23 / 4.0) - (
                              cgx * cgy * mg * vel[3] * cos26 / 4.0) - (cgx * cgy * mg * vel[4] * cos23 / 4.0) - (
                              cgx * cgy * mg * vel[4] * cos26 / 4.0) - (cgx * cgz * mg * vel[1] * cos05 / 2.0) + (
                              cgx * cgz * mg * vel[2] * cos05 / 2.0) - (cgx * cgz * mg * vel[3] * cos05 / 2.0) + (
                              cgx * cgz * mg * vel[4] * cos05 / 4.0) - (c3z * l3m * l2x * vel[1] * cos18) + (
                              c3z * l3m * l2x * vel[2] * cos18) - (c4x * l4m * l3z * vel[2] * cos12 / 2.0) + (
                              c4x * l4m * l3z * vel[3] * cos12 / 2.0) + (c5x * l5m * l4z * vel[3] * cos09 / 4.0) - (
                              c5x * l5m * l4z * vel[4] * cos09 / 4.0) + (cgx * l4z * mg * vel[3] * cos09 / 4.0) - (
                              cgx * l4z * mg * vel[4] * cos09 / 4.0) + (c5y * l5m * l5x * vel[1] * cos05 / 2.0) - (
                              c5y * l5m * l5x * vel[2] * cos05 / 2.0) + (c5y * l5m * l5x * vel[3] * cos05 / 2.0) - (
                              c5y * l5m * l5x * vel[4] * cos05 / 4.0) - (c5x * l5m * l5z * vel[1] * cos05 / 2.0) + (
                              c5x * l5m * l5z * vel[2] * cos05 / 2.0) - (c5x * l5m * l5z * vel[3] * cos05 / 2.0) + (
                              c5x * l5m * l5z * vel[4] * cos05 / 4.0) + (cgy * l5x * mg * vel[1] * cos05 / 2.0) - (
                              cgy * l5x * mg * vel[2] * cos05 / 2.0) + (cgy * l5x * mg * vel[3] * cos05 / 2.0) - (
                              cgy * l5x * mg * vel[4] * cos05 / 4.0) - (cgx * l5z * mg * vel[1] * cos05 / 2.0) + (
                              cgx * l5z * mg * vel[2] * cos05 / 2.0) - (cgx * l5z * mg * vel[3] * cos05 / 2.0) + (
                              cgx * l5z * mg * vel[4] * cos05 / 4.0) - (l4m * l2x * l4z * vel[1] * cos18) + (
                              l4m * l2x * l4z * vel[2] * cos18) - (l5m * l2x * l4z * vel[1] * cos18) + (
                              l5m * l2x * l4z * vel[2] * cos18) - (l5m * l5x * l3z * vel[2] * cos12 / 2.0) + (
                              l5m * l5x * l3z * vel[3] * cos12 / 2.0) - (l2x * l4z * mg * vel[1] * cos18) + (
                              l2x * l4z * mg * vel[2] * cos18) - (l5x * l3z * mg * vel[2] * cos12 / 2.0) + (
                              l5x * l3z * mg * vel[3] * cos12 / 2.0) - (c5y * c5z * l5m * vel[1] * sin12 / 2.0) + (
                              c5y * c5z * l5m * vel[2] * sin12 / 2.0) - (c5y * c5z * l5m * vel[3] * sin12 / 2.0) + (
                              c5y * c5z * l5m * vel[4] * sin12 / 4.0) - (cgy * cgz * mg * vel[1] * sin12 / 2.0) + (
                              cgy * cgz * mg * vel[2] * sin12 / 2.0) - (cgy * cgz * mg * vel[3] * sin12 / 2.0) + (
                              cgy * cgz * mg * vel[4] * sin12 / 4.0) - (c3x * l3m * l2x * vel[1] * sin30) + (
                              c3x * l3m * l2x * vel[2] * sin30) + (c5y * l5m * l4z * vel[3] * sin15 / 4.0) - (
                              c5y * l5m * l4z * vel[4] * sin15 / 4.0) - (c4z * l4m * l3z * vel[2] * sin11 / 2.0) + (
                              c4z * l4m * l3z * vel[3] * sin11 / 2.0) - (c5z * l5m * l3z * vel[2] * sin11 / 2.0) + (
                              c5z * l5m * l3z * vel[3] * sin11 / 2.0) + (cgy * l4z * mg * vel[3] * sin15 / 4.0) - (
                              cgy * l4z * mg * vel[4] * sin15 / 4.0) - (cgz * l3z * mg * vel[2] * sin11 / 2.0) + (
                              cgz * l3z * mg * vel[3] * sin11 / 2.0) - (c5x * l5m * l5x * vel[1] * sin12 / 2.0) + (
                              c5x * l5m * l5x * vel[2] * sin12 / 2.0) - (c5x * l5m * l5x * vel[3] * sin12 / 2.0) + (
                              c5x * l5m * l5x * vel[4] * sin12 / 4.0) - (c5y * l5m * l5z * vel[1] * sin12 / 2.0) + (
                              c5y * l5m * l5z * vel[2] * sin12 / 2.0) - (c5y * l5m * l5z * vel[3] * sin12 / 2.0) + (
                              c5y * l5m * l5z * vel[4] * sin12 / 4.0) - (cgx * l5x * mg * vel[1] * sin12 / 2.0) + (
                              cgx * l5x * mg * vel[2] * sin12 / 2.0) - (cgx * l5x * mg * vel[3] * sin12 / 2.0) + (
                              cgx * l5x * mg * vel[4] * sin12 / 4.0) - (cgy * l5z * mg * vel[1] * sin12 / 2.0) + (
                              cgy * l5z * mg * vel[2] * sin12 / 2.0) - (cgy * l5z * mg * vel[3] * sin12 / 2.0) + (
                              cgy * l5z * mg * vel[4] * sin12 / 4.0) - (l5m * l3z * l5z * vel[2] * sin11 / 2.0) + (
                              l5m * l3z * l5z * vel[3] * sin11 / 2.0) - (l3z * l5z * mg * vel[2] * sin11 / 2.0) + (
                              l3z * l5z * mg * vel[3] * sin11 / 2.0) + (c5x * c5y * l5m * vel[4] * cos11 / 2.0) - (
                              c2x * c2z * l2m * vel[1] * cos33) + (cgx * cgy * mg * vel[4] * cos11 / 2.0) - (
                              c4x * c4z * l4m * vel[1] * cos15) + (c4x * c4z * l4m * vel[2] * cos15) - (
                              c4x * c4z * l4m * vel[3] * cos15) - (c5z * l5m * l5x * vel[1] * cos15) + (
                              c5z * l5m * l5x * vel[2] * cos15) - (c5z * l5m * l5x * vel[3] * cos15) - (
                              cgz * l5x * mg * vel[1] * cos15) + (cgz * l5x * mg * vel[2] * cos15) - (
                              cgz * l5x * mg * vel[3] * cos15) - (l5m * l5x * l5z * vel[1] * cos15) + (
                              l5m * l5x * l5z * vel[2] * cos15) - (l5m * l5x * l5z * vel[3] * cos15) - (
                              l5x * l5z * mg * vel[1] * cos15) + (l5x * l5z * mg * vel[2] * cos15) - (
                              l5x * l5z * mg * vel[3] * cos15) + (c5z * l5m * l5z * vel[1] * sin05) - (
                              c5z * l5m * l5z * vel[2] * sin05) + (c5z * l5m * l5z * vel[3] * sin05) + (
                              cgz * l5z * mg * vel[1] * sin05) - (cgz * l5z * mg * vel[2] * sin05) + (
                              cgz * l5z * mg * vel[3] * sin05) + (c5y * l5m * l2x * vel[1] * cos00 / 2.0) - (
                              c5y * l5m * l2x * vel[2] * cos00 / 2.0) + (c5y * l5m * l2x * vel[3] * cos00 / 2.0) - (
                              c5y * l5m * l2x * vel[4] * cos00 / 2.0) - (c5x * l5m * l3z * vel[1] * cos25 / 2.0) + (
                              c5x * l5m * l3z * vel[2] * cos25 / 4.0) - (c5x * l5m * l3z * vel[3] * cos25 / 4.0) - (
                              c5x * l5m * l4z * vel[1] * cos24 / 2.0) - (c5x * l5m * l3z * vel[4] * cos25 / 4.0) + (
                              c5x * l5m * l4z * vel[2] * cos24 / 2.0) - (c5x * l5m * l4z * vel[3] * cos24 / 4.0) - (
                              c5x * l5m * l4z * vel[4] * cos24 / 4.0) + (cgy * l2x * mg * vel[1] * cos00 / 2.0) - (
                              cgy * l2x * mg * vel[2] * cos00 / 2.0) + (cgy * l2x * mg * vel[3] * cos00 / 2.0) - (
                              cgy * l2x * mg * vel[4] * cos00 / 2.0) - (cgx * l3z * mg * vel[1] * cos25 / 2.0) + (
                              cgx * l3z * mg * vel[2] * cos25 / 4.0) - (cgx * l3z * mg * vel[3] * cos25 / 4.0) - (
                              cgx * l4z * mg * vel[1] * cos24 / 2.0) - (cgx * l3z * mg * vel[4] * cos25 / 4.0) + (
                              cgx * l4z * mg * vel[2] * cos24 / 2.0) - (cgx * l4z * mg * vel[3] * cos24 / 4.0) - (
                              cgx * l4z * mg * vel[4] * cos24 / 4.0) - (c3x * c3z * l3m * vel[1] * cos31) + (
                              c3x * c3z * l3m * vel[2] * cos31) - (c5x * l5m * l2x * vel[1] * sin04 / 2.0) + (
                              c5x * l5m * l2x * vel[2] * sin04 / 2.0) - (c5x * l5m * l2x * vel[3] * sin04 / 2.0) + (
                              c5x * l5m * l2x * vel[4] * sin04 / 2.0) + (c5y * l5m * l3z * vel[1] * sin26 / 2.0) - (
                              c5y * l5m * l3z * vel[2] * sin26 / 4.0) + (c5y * l5m * l3z * vel[3] * sin26 / 4.0) + (
                              c5y * l5m * l4z * vel[1] * sin28 / 2.0) + (c5y * l5m * l3z * vel[4] * sin26 / 4.0) - (
                              c5y * l5m * l4z * vel[2] * sin28 / 2.0) + (c5y * l5m * l4z * vel[3] * sin28 / 4.0) + (
                              c5y * l5m * l4z * vel[4] * sin28 / 4.0) - (cgx * l2x * mg * vel[1] * sin04 / 2.0) + (
                              cgx * l2x * mg * vel[2] * sin04 / 2.0) - (cgx * l2x * mg * vel[3] * sin04 / 2.0) + (
                              cgx * l2x * mg * vel[4] * sin04 / 2.0) + (cgy * l3z * mg * vel[1] * sin26 / 2.0) - (
                              cgy * l3z * mg * vel[2] * sin26 / 4.0) + (cgy * l3z * mg * vel[3] * sin26 / 4.0) + (
                              cgy * l4z * mg * vel[1] * sin28 / 2.0) + (cgy * l3z * mg * vel[4] * sin26 / 4.0) - (
                              cgy * l4z * mg * vel[2] * sin28 / 2.0) + (cgy * l4z * mg * vel[3] * sin28 / 4.0) + (
                              cgy * l4z * mg * vel[4] * sin28 / 4.0) - (c3x * l3m * l3z * vel[1] * cos32) + (
                              c3x * l3m * l3z * vel[2] * cos32 / 2.0) + (c3z * l3m * l3z * vel[1] * sin29) - (
                              c3z * l3m * l3z * vel[2] * sin29 / 2.0) - (c4z * l4m * l2x * vel[1] * cos04) + (
                              c4z * l4m * l2x * vel[2] * cos04) - (c4z * l4m * l2x * vel[3] * cos04) - (
                              c5z * l5m * l2x * vel[1] * cos04) + (c5z * l5m * l2x * vel[2] * cos04) - (
                              c5z * l5m * l2x * vel[3] * cos04) - (c5x * l5m * l3z * vel[2] * cos14 / 4.0) + (
                              c5x * l5m * l3z * vel[3] * cos14 / 4.0) - (c5x * l5m * l3z * vel[4] * cos14 / 4.0) - (
                              cgz * l2x * mg * vel[1] * cos04) + (cgz * l2x * mg * vel[2] * cos04) - (
                              cgz * l2x * mg * vel[3] * cos04) - (cgx * l3z * mg * vel[2] * cos14 / 4.0) + (
                              cgx * l3z * mg * vel[3] * cos14 / 4.0) - (cgx * l3z * mg * vel[4] * cos14 / 4.0) + (
                              l4m * l3z * l4z * vel[1] * sin29) - (l4m * l3z * l4z * vel[2] * sin29 / 2.0) + (
                              l5m * l3z * l4z * vel[1] * sin29) - (l5m * l3z * l4z * vel[2] * sin29 / 2.0) + (
                              l3z * l4z * mg * vel[1] * sin29) - (l3z * l4z * mg * vel[2] * sin29 / 2.0) - (
                              l5m * l2x * l5z * vel[1] * cos04) + (l5m * l2x * l5z * vel[2] * cos04) - (
                              l5m * l2x * l5z * vel[3] * cos04) - (l2x * l5z * mg * vel[1] * cos04) + (
                              l2x * l5z * mg * vel[2] * cos04) - (l2x * l5z * mg * vel[3] * cos04) - (
                              c4x * l4m * l2x * vel[1] * sin06) + (c4x * l4m * l2x * vel[2] * sin06) - (
                              c4x * l4m * l2x * vel[3] * sin06) + (c5y * l5m * l3z * vel[2] * sin16 / 4.0) - (
                              c5y * l5m * l3z * vel[3] * sin16 / 4.0) + (c5y * l5m * l3z * vel[4] * sin16 / 4.0) + (
                              cgy * l3z * mg * vel[2] * sin16 / 4.0) - (cgy * l3z * mg * vel[3] * sin16 / 4.0) + (
                              cgy * l3z * mg * vel[4] * sin16 / 4.0) - (l5m * l2x * l5x * vel[1] * sin06) + (
                              l5m * l2x * l5x * vel[2] * sin06) - (l5m * l2x * l5x * vel[3] * sin06) - (
                              l2x * l5x * mg * vel[1] * sin06) + (l2x * l5x * mg * vel[2] * sin06) - (
                              l2x * l5x * mg * vel[3] * sin06) - (c5x * c5z * l5m * vel[1] * cos06 / 2.0) + (
                              c5x * c5z * l5m * vel[2] * cos06 / 2.0) - (c5x * c5z * l5m * vel[3] * cos06 / 2.0) - (
                              c5x * c5z * l5m * vel[4] * cos06 / 4.0) - (cgx * cgz * mg * vel[1] * cos06 / 2.0) + (
                              cgx * cgz * mg * vel[2] * cos06 / 2.0) - (cgx * cgz * mg * vel[3] * cos06 / 2.0) - (
                              cgx * cgz * mg * vel[4] * cos06 / 4.0) + (c5x * l5m * l4z * vel[3] * cos07 / 4.0) + (
                              c5x * l5m * l4z * vel[4] * cos07 / 4.0) + (cgx * l4z * mg * vel[3] * cos07 / 4.0) + (
                              cgx * l4z * mg * vel[4] * cos07 / 4.0) - (c5y * l5m * l5x * vel[1] * cos06 / 2.0) + (
                              c5y * l5m * l5x * vel[2] * cos06 / 2.0) - (c5y * l5m * l5x * vel[3] * cos06 / 2.0) - (
                              c5y * l5m * l5x * vel[4] * cos06 / 4.0) - (c5x * l5m * l3z * vel[1] * cos22 / 2.0) + (
                              c5x * l5m * l3z * vel[2] * cos22 / 4.0) - (c5x * l5m * l3z * vel[3] * cos22 / 4.0) - (
                              c5x * l5m * l4z * vel[1] * cos21 / 2.0) + (c5x * l5m * l3z * vel[4] * cos22 / 4.0) + (
                              c5x * l5m * l4z * vel[2] * cos21 / 2.0) - (c5x * l5m * l4z * vel[3] * cos21 / 4.0) - (
                              c5x * l5m * l5z * vel[1] * cos06 / 2.0) + (c5x * l5m * l4z * vel[4] * cos21 / 4.0) + (
                              c5x * l5m * l5z * vel[2] * cos06 / 2.0) - (c5x * l5m * l5z * vel[3] * cos06 / 2.0) - (
                              c5x * l5m * l5z * vel[4] * cos06 / 4.0) - (cgy * l5x * mg * vel[1] * cos06 / 2.0) + (
                              cgy * l5x * mg * vel[2] * cos06 / 2.0) - (cgy * l5x * mg * vel[3] * cos06 / 2.0) - (
                              cgy * l5x * mg * vel[4] * cos06 / 4.0) - (cgx * l3z * mg * vel[1] * cos22 / 2.0) + (
                              cgx * l3z * mg * vel[2] * cos22 / 4.0) - (cgx * l3z * mg * vel[3] * cos22 / 4.0) - (
                              cgx * l4z * mg * vel[1] * cos21 / 2.0) + (cgx * l3z * mg * vel[4] * cos22 / 4.0) + (
                              cgx * l4z * mg * vel[2] * cos21 / 2.0) - (cgx * l4z * mg * vel[3] * cos21 / 4.0) - (
                              cgx * l5z * mg * vel[1] * cos06 / 2.0) + (cgx * l4z * mg * vel[4] * cos21 / 4.0) + (
                              cgx * l5z * mg * vel[2] * cos06 / 2.0) - (cgx * l5z * mg * vel[3] * cos06 / 2.0) - (
                              cgx * l5z * mg * vel[4] * cos06 / 4.0) + (c5y * c5z * l5m * vel[1] * sin13 / 2.0) - (
                              c5y * c5z * l5m * vel[2] * sin13 / 2.0) + (c5y * c5z * l5m * vel[3] * sin13 / 2.0) + (
                              c5y * c5z * l5m * vel[4] * sin13 / 4.0) + (cgy * cgz * mg * vel[1] * sin13 / 2.0) - (
                              cgy * cgz * mg * vel[2] * sin13 / 2.0) + (cgy * cgz * mg * vel[3] * sin13 / 2.0) + (
                              cgy * cgz * mg * vel[4] * sin13 / 4.0) - (c5y * l5m * l4z * vel[3] * sin14 / 4.0) - (
                              c5y * l5m * l4z * vel[4] * sin14 / 4.0) - (cgy * l4z * mg * vel[3] * sin14 / 4.0) - (
                              cgy * l4z * mg * vel[4] * sin14 / 4.0) - (c5x * l5m * l5x * vel[1] * sin13 / 2.0) + (
                              c5x * l5m * l5x * vel[2] * sin13 / 2.0) - (c5x * l5m * l5x * vel[3] * sin13 / 2.0) - (
                              c5x * l5m * l5x * vel[4] * sin13 / 4.0) - (c5y * l5m * l3z * vel[1] * sin27 / 2.0) + (
                              c5y * l5m * l3z * vel[2] * sin27 / 4.0) - (c5y * l5m * l3z * vel[3] * sin27 / 4.0) - (
                              c5y * l5m * l4z * vel[1] * sin25 / 2.0) + (c5y * l5m * l3z * vel[4] * sin27 / 4.0) + (
                              c5y * l5m * l4z * vel[2] * sin25 / 2.0) - (c5y * l5m * l4z * vel[3] * sin25 / 4.0) + (
                              c5y * l5m * l5z * vel[1] * sin13 / 2.0) + (c5y * l5m * l4z * vel[4] * sin25 / 4.0) - (
                              c5y * l5m * l5z * vel[2] * sin13 / 2.0) + (c5y * l5m * l5z * vel[3] * sin13 / 2.0) + (
                              c5y * l5m * l5z * vel[4] * sin13 / 4.0) - (cgx * l5x * mg * vel[1] * sin13 / 2.0) + (
                              cgx * l5x * mg * vel[2] * sin13 / 2.0) - (cgx * l5x * mg * vel[3] * sin13 / 2.0) - (
                              cgx * l5x * mg * vel[4] * sin13 / 4.0) - (cgy * l3z * mg * vel[1] * sin27 / 2.0) + (
                              cgy * l3z * mg * vel[2] * sin27 / 4.0) - (cgy * l3z * mg * vel[3] * sin27 / 4.0) - (
                              cgy * l4z * mg * vel[1] * sin25 / 2.0) + (cgy * l3z * mg * vel[4] * sin27 / 4.0) + (
                              cgy * l4z * mg * vel[2] * sin25 / 2.0) - (cgy * l4z * mg * vel[3] * sin25 / 4.0) + (
                              cgy * l5z * mg * vel[1] * sin13 / 2.0) + (cgy * l4z * mg * vel[4] * sin25 / 4.0) - (
                              cgy * l5z * mg * vel[2] * sin13 / 2.0) + (cgy * l5z * mg * vel[3] * sin13 / 2.0) + (
                              cgy * l5z * mg * vel[4] * sin13 / 4.0) - (c5y * l5m * l5x * vel[4] * cos02 / 2.0) - (
                              c2z * l2m * l2x * vel[1] * cos27) - (c3x * l3m * l3z * vel[2] * cos30 / 2.0) + (
                              c4x * l4m * l4z * vel[3] * cos03 / 2.0) - (cgy * l5x * mg * vel[4] * cos02 / 2.0) - (
                              c4x * l4m * l3z * vel[1] * cos20) + (c4x * l4m * l3z * vel[2] * cos20 / 2.0) - (
                              c4x * l4m * l3z * vel[3] * cos20 / 2.0) - (c4x * l4m * l4z * vel[1] * cos19) + (
                              c4x * l4m * l4z * vel[2] * cos19) - (c4x * l4m * l4z * vel[3] * cos19 / 2.0) - (
                              c5x * l5m * l3z * vel[2] * cos13 / 4.0) + (c5x * l5m * l3z * vel[3] * cos13 / 4.0) + (
                              c5x * l5m * l3z * vel[4] * cos13 / 4.0) - (cgx * l3z * mg * vel[2] * cos13 / 4.0) + (
                              cgx * l3z * mg * vel[3] * cos13 / 4.0) + (cgx * l3z * mg * vel[4] * cos13 / 4.0) - (
                              l3m * l2x * l3z * vel[1] * cos27) - (l4m * l2x * l3z * vel[1] * cos27) - (
                              l5m * l2x * l3z * vel[1] * cos27) + (l5m * l5x * l4z * vel[3] * cos03 / 2.0) - (
                              l2x * l3z * mg * vel[1] * cos27) + (l5x * l4z * mg * vel[3] * cos03 / 2.0) - (
                              l5m * l5x * l3z * vel[1] * cos20) + (l5m * l5x * l3z * vel[2] * cos20 / 2.0) - (
                              l5m * l5x * l3z * vel[3] * cos20 / 2.0) - (l5m * l5x * l4z * vel[1] * cos19) + (
                              l5m * l5x * l4z * vel[2] * cos19) - (l5m * l5x * l4z * vel[3] * cos19 / 2.0) - (
                              c2x * l2m * l2x * vel[1] * sin33) - (c5x * l5m * l5x * vel[4] * sin07 / 2.0) - (
                              c3z * l3m * l3z * vel[2] * sin24 / 2.0) - (c4z * l4m * l4z * vel[3] * sin08 / 2.0) - (
                              c5z * l5m * l4z * vel[3] * sin08 / 2.0) - (l5x * l3z * mg * vel[1] * cos20) + (
                              l5x * l3z * mg * vel[2] * cos20 / 2.0) - (l5x * l3z * mg * vel[3] * cos20 / 2.0) - (
                              l5x * l4z * mg * vel[1] * cos19) + (l5x * l4z * mg * vel[2] * cos19) - (
                              l5x * l4z * mg * vel[3] * cos19 / 2.0) - (cgx * l5x * mg * vel[4] * sin07 / 2.0) - (
                              cgz * l4z * mg * vel[3] * sin08 / 2.0) - (c5y * l5m * l3z * vel[2] * sin17 / 4.0) + (
                              c5y * l5m * l3z * vel[3] * sin17 / 4.0) + (c5y * l5m * l3z * vel[4] * sin17 / 4.0) + (
                              c4z * l4m * l3z * vel[1] * sin19) - (c4z * l4m * l3z * vel[2] * sin19 / 2.0) + (
                              c4z * l4m * l3z * vel[3] * sin19 / 2.0) + (c4z * l4m * l4z * vel[1] * sin20) + (
                              c5z * l5m * l3z * vel[1] * sin19) - (c4z * l4m * l4z * vel[2] * sin20) - (
                              c5z * l5m * l3z * vel[2] * sin19 / 2.0) + (c4z * l4m * l4z * vel[3] * sin20 / 2.0) + (
                              c5z * l5m * l3z * vel[3] * sin19 / 2.0) + (c5z * l5m * l4z * vel[1] * sin20) - (
                              c5z * l5m * l4z * vel[2] * sin20) + (c5z * l5m * l4z * vel[3] * sin20 / 2.0) - (
                              cgy * l3z * mg * vel[2] * sin17 / 4.0) + (cgy * l3z * mg * vel[3] * sin17 / 4.0) + (
                              cgy * l3z * mg * vel[4] * sin17 / 4.0) + (cgz * l3z * mg * vel[1] * sin19) - (
                              cgz * l3z * mg * vel[2] * sin19 / 2.0) + (cgz * l3z * mg * vel[3] * sin19 / 2.0) + (
                              cgz * l4z * mg * vel[1] * sin20) - (cgz * l4z * mg * vel[2] * sin20) + (
                              cgz * l4z * mg * vel[3] * sin20 / 2.0) - (l4m * l3z * l4z * vel[2] * sin24 / 2.0) - (
                              l5m * l3z * l4z * vel[2] * sin24 / 2.0) - (l5m * l4z * l5z * vel[3] * sin08 / 2.0)
        C[0, 1] = (I5x * vel[0] * sin10 / 8.0) + (I5x * vel[0] * sin09 / 8.0) - (I5y * vel[0] * sin10 / 8.0) - (
                    I5y * vel[0] * sin09 / 8.0) + (Igx * vel[0] * sin10 / 8.0) + (Igx * vel[0] * sin09 / 8.0) - (
                              Igy * vel[0] * sin10 / 8.0) - (Igy * vel[0] * sin09 / 8.0) + (
                              I2x * vel[0] * sin23 / 2.0) - (I2z * vel[0] * sin23 / 2.0) + (
                              I4x * vel[0] * sin05 / 2.0) + (I5x * vel[0] * sin05 / 4.0) + (
                              I5y * vel[0] * sin05 / 4.0) - (I4z * vel[0] * sin05 / 2.0) - (
                              I5z * vel[0] * sin05 / 2.0) + (Igx * vel[0] * sin05 / 4.0) + (
                              Igy * vel[0] * sin05 / 4.0) - (Igz * vel[0] * sin05 / 2.0) + (
                              I5x * vel[1] * sin01 / 4.0) - (I5x * vel[1] * sin02 / 4.0) - (
                              I5x * vel[2] * sin01 / 4.0) + (I5x * vel[2] * sin02 / 4.0) + (
                              I5x * vel[3] * sin01 / 4.0) - (I5x * vel[3] * sin02 / 4.0) - (
                              I5x * vel[4] * sin01 / 4.0) - (I5x * vel[4] * sin02 / 4.0) - (
                              I5y * vel[1] * sin01 / 4.0) + (I5y * vel[1] * sin02 / 4.0) + (
                              I5y * vel[2] * sin01 / 4.0) - (I5y * vel[2] * sin02 / 4.0) - (
                              I5y * vel[3] * sin01 / 4.0) + (I5y * vel[3] * sin02 / 4.0) + (
                              I5y * vel[4] * sin01 / 4.0) + (I5y * vel[4] * sin02 / 4.0) + (
                              Igx * vel[1] * sin01 / 4.0) - (Igx * vel[1] * sin02 / 4.0) - (
                              Igx * vel[2] * sin01 / 4.0) + (Igx * vel[2] * sin02 / 4.0) + (
                              Igx * vel[3] * sin01 / 4.0) - (Igx * vel[3] * sin02 / 4.0) - (
                              Igx * vel[4] * sin01 / 4.0) - (Igx * vel[4] * sin02 / 4.0) - (
                              Igy * vel[1] * sin01 / 4.0) + (Igy * vel[1] * sin02 / 4.0) + (
                              Igy * vel[2] * sin01 / 4.0) - (Igy * vel[2] * sin02 / 4.0) - (
                              Igy * vel[3] * sin01 / 4.0) + (Igy * vel[3] * sin02 / 4.0) + (
                              Igy * vel[4] * sin01 / 4.0) + (Igy * vel[4] * sin02 / 4.0) + (
                              I3x * vel[0] * sin18 / 2.0) - (I3z * vel[0] * sin18 / 2.0) - (
                              I5z * vel[4] * sin06 / 2.0) - (Igz * vel[4] * sin06 / 2.0) - (
                              pow(c5x, 2.0) * l5m * vel[0] * sin10 / 8.0) - (
                              pow(c5x, 2.0) * l5m * vel[0] * sin09 / 8.0) + (
                              pow(c5y, 2.0) * l5m * vel[0] * sin10 / 8.0) + (
                              pow(c5y, 2.0) * l5m * vel[0] * sin09 / 8.0) - (
                              pow(cgx, 2.0) * mg * vel[0] * sin10 / 8.0) - (
                              pow(cgx, 2.0) * mg * vel[0] * sin09 / 8.0) + (
                              pow(cgy, 2.0) * mg * vel[0] * sin10 / 8.0) + (
                              pow(cgy, 2.0) * mg * vel[0] * sin09 / 8.0) - (
                              pow(c2x, 2.0) * l2m * vel[0] * sin23 / 2.0) + (
                              pow(c2z, 2.0) * l2m * vel[0] * sin23 / 2.0) - (
                              pow(c4x, 2.0) * l4m * vel[0] * sin05 / 2.0) - (
                              pow(c5x, 2.0) * l5m * vel[0] * sin05 / 4.0) - (
                              pow(c5y, 2.0) * l5m * vel[0] * sin05 / 4.0) + (
                              pow(c4z, 2.0) * l4m * vel[0] * sin05 / 2.0) + (
                              pow(c5z, 2.0) * l5m * vel[0] * sin05 / 2.0) - (
                              pow(cgx, 2.0) * mg * vel[0] * sin05 / 4.0) - (
                              pow(cgy, 2.0) * mg * vel[0] * sin05 / 4.0) + (
                              pow(cgz, 2.0) * mg * vel[0] * sin05 / 2.0) + (
                              l3m * pow(l3z, 2.0) * vel[0] * sin23 / 2.0) + (
                              l4m * pow(l3z, 2.0) * vel[0] * sin23 / 2.0) + (
                              l5m * pow(l3z, 2.0) * vel[0] * sin23 / 2.0) + (
                              pow(l3z, 2.0) * mg * vel[0] * sin23 / 2.0) - (
                              l5m * pow(l5x, 2.0) * vel[0] * sin05 / 2.0) + (
                              l5m * pow(l5z, 2.0) * vel[0] * sin05 / 2.0) - (
                              pow(l5x, 2.0) * mg * vel[0] * sin05 / 2.0) + (
                              pow(l5z, 2.0) * mg * vel[0] * sin05 / 2.0) - (
                              pow(c5x, 2.0) * l5m * vel[1] * sin01 / 4.0) + (
                              pow(c5x, 2.0) * l5m * vel[1] * sin02 / 4.0) + (
                              pow(c5x, 2.0) * l5m * vel[2] * sin01 / 4.0) - (
                              pow(c5x, 2.0) * l5m * vel[2] * sin02 / 4.0) - (
                              pow(c5x, 2.0) * l5m * vel[3] * sin01 / 4.0) + (
                              pow(c5x, 2.0) * l5m * vel[3] * sin02 / 4.0) + (
                              pow(c5x, 2.0) * l5m * vel[4] * sin01 / 4.0) + (
                              pow(c5x, 2.0) * l5m * vel[4] * sin02 / 4.0) + (
                              pow(c5y, 2.0) * l5m * vel[1] * sin01 / 4.0) - (
                              pow(c5y, 2.0) * l5m * vel[1] * sin02 / 4.0) - (
                              pow(c5y, 2.0) * l5m * vel[2] * sin01 / 4.0) + (
                              pow(c5y, 2.0) * l5m * vel[2] * sin02 / 4.0) + (
                              pow(c5y, 2.0) * l5m * vel[3] * sin01 / 4.0) - (
                              pow(c5y, 2.0) * l5m * vel[3] * sin02 / 4.0) - (
                              pow(c5y, 2.0) * l5m * vel[4] * sin01 / 4.0) - (
                              pow(c5y, 2.0) * l5m * vel[4] * sin02 / 4.0) - (
                              pow(cgx, 2.0) * mg * vel[1] * sin01 / 4.0) + (
                              pow(cgx, 2.0) * mg * vel[1] * sin02 / 4.0) + (
                              pow(cgx, 2.0) * mg * vel[2] * sin01 / 4.0) - (
                              pow(cgx, 2.0) * mg * vel[2] * sin02 / 4.0) - (
                              pow(cgx, 2.0) * mg * vel[3] * sin01 / 4.0) + (
                              pow(cgx, 2.0) * mg * vel[3] * sin02 / 4.0) + (
                              pow(cgx, 2.0) * mg * vel[4] * sin01 / 4.0) + (
                              pow(cgx, 2.0) * mg * vel[4] * sin02 / 4.0) + (
                              pow(cgy, 2.0) * mg * vel[1] * sin01 / 4.0) - (
                              pow(cgy, 2.0) * mg * vel[1] * sin02 / 4.0) - (
                              pow(cgy, 2.0) * mg * vel[2] * sin01 / 4.0) + (
                              pow(cgy, 2.0) * mg * vel[2] * sin02 / 4.0) + (
                              pow(cgy, 2.0) * mg * vel[3] * sin01 / 4.0) - (
                              pow(cgy, 2.0) * mg * vel[3] * sin02 / 4.0) - (
                              pow(cgy, 2.0) * mg * vel[4] * sin01 / 4.0) - (
                              pow(cgy, 2.0) * mg * vel[4] * sin02 / 4.0) - (
                              pow(c3x, 2.0) * l3m * vel[0] * sin18 / 2.0) + (
                              pow(c3z, 2.0) * l3m * vel[0] * sin18 / 2.0) + (
                              l4m * pow(l4z, 2.0) * vel[0] * sin18 / 2.0) + (
                              l5m * pow(l4z, 2.0) * vel[0] * sin18 / 2.0) + (
                              pow(l4z, 2.0) * mg * vel[0] * sin18 / 2.0) - (
                              pow(c5x, 2.0) * l5m * vel[4] * sin06 / 2.0) - (
                              pow(c5y, 2.0) * l5m * vel[4] * sin06 / 2.0) - (
                              pow(cgx, 2.0) * mg * vel[4] * sin06 / 2.0) - (
                              pow(cgy, 2.0) * mg * vel[4] * sin06 / 2.0) - (c5y * l5m * l2x * vel[0] * cos01 / 2.0) + (
                              c5y * l5m * l5x * vel[1] * cos01 / 2.0) - (c5y * l5m * l5x * vel[2] * cos01 / 2.0) + (
                              c5y * l5m * l5x * vel[3] * cos01 / 2.0) + (c5x * l5m * l5z * vel[1] * cos01 / 2.0) - (
                              c5x * l5m * l5z * vel[2] * cos01 / 2.0) + (c5x * l5m * l5z * vel[3] * cos01 / 2.0) - (
                              cgy * l2x * mg * vel[0] * cos01 / 2.0) + (cgy * l5x * mg * vel[1] * cos01 / 2.0) - (
                              cgy * l5x * mg * vel[2] * cos01 / 2.0) + (cgy * l5x * mg * vel[3] * cos01 / 2.0) + (
                              cgx * l5z * mg * vel[1] * cos01 / 2.0) - (cgx * l5z * mg * vel[2] * cos01 / 2.0) + (
                              cgx * l5z * mg * vel[3] * cos01 / 2.0) + (l5m * l3z * l5z * vel[0] * sin19) + (
                              l5m * l4z * l5z * vel[0] * sin20) + (l3z * l5z * mg * vel[0] * sin19) + (
                              l4z * l5z * mg * vel[0] * sin20) - (c5y * c5z * l5m * vel[1] * sin03 / 2.0) + (
                              c5y * c5z * l5m * vel[2] * sin03 / 2.0) - (c5y * c5z * l5m * vel[3] * sin03 / 2.0) - (
                              cgy * cgz * mg * vel[1] * sin03 / 2.0) + (cgy * cgz * mg * vel[2] * sin03 / 2.0) - (
                              cgy * cgz * mg * vel[3] * sin03 / 2.0) + (c3x * c3y * l3m * vel[1] * cos18) - (
                              c3x * c3y * l3m * vel[2] * cos18) - (c5x * l5m * l2x * vel[0] * sin03 / 2.0) + (
                              c5x * l5m * l5x * vel[1] * sin03 / 2.0) - (c5x * l5m * l5x * vel[2] * sin03 / 2.0) + (
                              c5x * l5m * l5x * vel[3] * sin03 / 2.0) - (c5y * l5m * l5z * vel[1] * sin03 / 2.0) + (
                              c5y * l5m * l5z * vel[2] * sin03 / 2.0) - (c5y * l5m * l5z * vel[3] * sin03 / 2.0) - (
                              cgx * l2x * mg * vel[0] * sin03 / 2.0) + (cgx * l5x * mg * vel[1] * sin03 / 2.0) - (
                              cgx * l5x * mg * vel[2] * sin03 / 2.0) + (cgx * l5x * mg * vel[3] * sin03 / 2.0) - (
                              cgy * l5z * mg * vel[1] * sin03 / 2.0) + (cgy * l5z * mg * vel[2] * sin03 / 2.0) - (
                              cgy * l5z * mg * vel[3] * sin03 / 2.0) + (c5x * c5y * l5m * vel[0] * cos23 / 4.0) - (
                              c5x * c5y * l5m * vel[0] * cos26 / 4.0) - (c5x * c5z * l5m * vel[0] * cos05 / 2.0) + (
                              cgx * cgy * mg * vel[0] * cos23 / 4.0) - (cgx * cgy * mg * vel[0] * cos26 / 4.0) - (
                              cgx * cgz * mg * vel[0] * cos05 / 2.0) - (c3z * l3m * l2x * vel[0] * cos18) - (
                              c5x * l5m * l3z * vel[1] * cos29 / 2.0) - (cgx * l3z * mg * vel[1] * cos29 / 2.0) + (
                              c5y * l5m * l5x * vel[0] * cos05 / 2.0) - (c5x * l5m * l5z * vel[0] * cos05 / 2.0) - (
                              c3y * c3z * l3m * vel[1] * sin30) + (c3y * c3z * l3m * vel[2] * sin30) + (
                              cgy * l5x * mg * vel[0] * cos05 / 2.0) - (cgx * l5z * mg * vel[0] * cos05 / 2.0) - (
                              l4m * l2x * l4z * vel[0] * cos18) - (l5m * l2x * l4z * vel[0] * cos18) - (
                              l2x * l4z * mg * vel[0] * cos18) - (c5y * c5z * l5m * vel[0] * sin12 / 2.0) - (
                              cgy * cgz * mg * vel[0] * sin12 / 2.0) - (c3x * l3m * l2x * vel[0] * sin30) - (
                              c4y * l4m * l4z * vel[1] * sin30) + (c4y * l4m * l4z * vel[2] * sin30) - (
                              c5y * l5m * l3z * vel[1] * sin32 / 2.0) - (cgy * l3z * mg * vel[1] * sin32 / 2.0) - (
                              c5x * l5m * l5x * vel[0] * sin12 / 2.0) - (c5y * l5m * l5z * vel[0] * sin12 / 2.0) - (
                              cgx * l5x * mg * vel[0] * sin12 / 2.0) - (cgy * l5z * mg * vel[0] * sin12 / 2.0) - (
                              c2x * c2z * l2m * vel[0] * cos33) - (c4x * c4z * l4m * vel[0] * cos15) - (
                              c5z * l5m * l5x * vel[0] * cos15) - (cgz * l5x * mg * vel[0] * cos15) - (
                              l5m * l5x * l5z * vel[0] * cos15) - (l5x * l5z * mg * vel[0] * cos15) + (
                              c5x * c5y * l5m * vel[1] * cos08 / 2.0) + (c5x * c5y * l5m * vel[1] * cos10 / 2.0) - (
                              c5x * c5y * l5m * vel[2] * cos08 / 2.0) - (c5x * c5y * l5m * vel[2] * cos10 / 2.0) + (
                              c5x * c5y * l5m * vel[3] * cos08 / 2.0) + (c5x * c5y * l5m * vel[3] * cos10 / 2.0) - (
                              c5x * c5y * l5m * vel[4] * cos08 / 2.0) + (c5x * c5y * l5m * vel[4] * cos10 / 2.0) - (
                              c5x * c5z * l5m * vel[1] * cos00 / 2.0) + (c5x * c5z * l5m * vel[2] * cos00 / 2.0) - (
                              c5x * c5z * l5m * vel[3] * cos00 / 2.0) + (cgx * cgy * mg * vel[1] * cos08 / 2.0) + (
                              cgx * cgy * mg * vel[1] * cos10 / 2.0) - (cgx * cgy * mg * vel[2] * cos08 / 2.0) - (
                              cgx * cgy * mg * vel[2] * cos10 / 2.0) + (cgx * cgy * mg * vel[3] * cos08 / 2.0) + (
                              cgx * cgy * mg * vel[3] * cos10 / 2.0) - (cgx * cgy * mg * vel[4] * cos08 / 2.0) + (
                              cgx * cgy * mg * vel[4] * cos10 / 2.0) - (cgx * cgz * mg * vel[1] * cos00 / 2.0) + (
                              cgx * cgz * mg * vel[2] * cos00 / 2.0) - (cgx * cgz * mg * vel[3] * cos00 / 2.0) + (
                              c5z * l5m * l5z * vel[0] * sin05) + (cgz * l5z * mg * vel[0] * sin05) + (
                              c5y * l5m * l2x * vel[0] * cos00 / 2.0) + (c5y * l5m * l5x * vel[1] * cos00 / 2.0) - (
                              c5y * l5m * l5x * vel[2] * cos00 / 2.0) + (c5y * l5m * l5x * vel[3] * cos00 / 2.0) - (
                              c5x * l5m * l3z * vel[0] * cos25 / 2.0) - (c5x * l5m * l4z * vel[0] * cos24 / 2.0) - (
                              c5x * l5m * l5z * vel[1] * cos00 / 2.0) + (c5x * l5m * l5z * vel[2] * cos00 / 2.0) - (
                              c5x * l5m * l5z * vel[3] * cos00 / 2.0) + (cgy * l2x * mg * vel[0] * cos00 / 2.0) + (
                              cgy * l5x * mg * vel[1] * cos00 / 2.0) - (cgy * l5x * mg * vel[2] * cos00 / 2.0) + (
                              cgy * l5x * mg * vel[3] * cos00 / 2.0) - (cgx * l3z * mg * vel[0] * cos25 / 2.0) - (
                              cgx * l4z * mg * vel[0] * cos24 / 2.0) - (cgx * l5z * mg * vel[1] * cos00 / 2.0) + (
                              cgx * l5z * mg * vel[2] * cos00 / 2.0) - (cgx * l5z * mg * vel[3] * cos00 / 2.0) - (
                              c5y * c5z * l5m * vel[1] * sin04 / 2.0) + (c5y * c5z * l5m * vel[2] * sin04 / 2.0) - (
                              c5y * c5z * l5m * vel[3] * sin04 / 2.0) - (cgy * cgz * mg * vel[1] * sin04 / 2.0) + (
                              cgy * cgz * mg * vel[2] * sin04 / 2.0) - (cgy * cgz * mg * vel[3] * sin04 / 2.0) - (
                              c3x * c3z * l3m * vel[0] * cos31) - (c5x * l5m * l2x * vel[0] * sin04 / 2.0) - (
                              c5x * l5m * l5x * vel[1] * sin04 / 2.0) + (c5x * l5m * l5x * vel[2] * sin04 / 2.0) - (
                              c5x * l5m * l5x * vel[3] * sin04 / 2.0) + (c5y * l5m * l3z * vel[0] * sin26 / 2.0) + (
                              c5y * l5m * l4z * vel[0] * sin28 / 2.0) - (c5y * l5m * l5z * vel[1] * sin04 / 2.0) + (
                              c5y * l5m * l5z * vel[2] * sin04 / 2.0) - (c5y * l5m * l5z * vel[3] * sin04 / 2.0) - (
                              cgx * l2x * mg * vel[0] * sin04 / 2.0) - (cgx * l5x * mg * vel[1] * sin04 / 2.0) + (
                              cgx * l5x * mg * vel[2] * sin04 / 2.0) - (cgx * l5x * mg * vel[3] * sin04 / 2.0) + (
                              cgy * l3z * mg * vel[0] * sin26 / 2.0) + (cgy * l4z * mg * vel[0] * sin28 / 2.0) - (
                              cgy * l5z * mg * vel[1] * sin04 / 2.0) + (cgy * l5z * mg * vel[2] * sin04 / 2.0) - (
                              cgy * l5z * mg * vel[3] * sin04 / 2.0) - (c3x * l3m * l3z * vel[0] * cos32) + (
                              c4x * c4y * l4m * vel[1] * cos04) - (c4x * c4y * l4m * vel[2] * cos04) + (
                              c4x * c4y * l4m * vel[3] * cos04) + (c3z * l3m * l3z * vel[0] * sin29) - (
                              c4z * l4m * l2x * vel[0] * cos04) - (c5z * l5m * l2x * vel[0] * cos04) + (
                              c5x * l5m * l4z * vel[1] * cos16 / 2.0) - (c5x * l5m * l4z * vel[2] * cos16 / 2.0) - (
                              cgz * l2x * mg * vel[0] * cos04) + (cgx * l4z * mg * vel[1] * cos16 / 2.0) - (
                              cgx * l4z * mg * vel[2] * cos16 / 2.0) + (l4m * l3z * l4z * vel[0] * sin29) + (
                              l5m * l3z * l4z * vel[0] * sin29) + (l3z * l4z * mg * vel[0] * sin29) - (
                              c4y * c4z * l4m * vel[1] * sin06) + (c4y * c4z * l4m * vel[2] * sin06) - (
                              c4y * c4z * l4m * vel[3] * sin06) - (l5m * l2x * l5z * vel[0] * cos04) - (
                              l2x * l5z * mg * vel[0] * cos04) - (c4x * l4m * l2x * vel[0] * sin06) - (
                              c5y * l5m * l4z * vel[1] * sin21 / 2.0) + (c5y * l5m * l4z * vel[2] * sin21 / 2.0) - (
                              cgy * l4z * mg * vel[1] * sin21 / 2.0) + (cgy * l4z * mg * vel[2] * sin21 / 2.0) - (
                              l5m * l2x * l5x * vel[0] * sin06) - (l2x * l5x * mg * vel[0] * sin06) - (
                              c5x * c5z * l5m * vel[0] * cos06 / 2.0) - (cgx * cgz * mg * vel[0] * cos06 / 2.0) + (
                              c5x * l5m * l3z * vel[1] * cos28 / 2.0) + (cgx * l3z * mg * vel[1] * cos28 / 2.0) - (
                              c5y * l5m * l5x * vel[0] * cos06 / 2.0) - (c5x * l5m * l3z * vel[0] * cos22 / 2.0) - (
                              c5x * l5m * l4z * vel[0] * cos21 / 2.0) - (c5x * l5m * l5z * vel[0] * cos06 / 2.0) - (
                              cgy * l5x * mg * vel[0] * cos06 / 2.0) - (cgx * l3z * mg * vel[0] * cos22 / 2.0) - (
                              cgx * l4z * mg * vel[0] * cos21 / 2.0) - (cgx * l5z * mg * vel[0] * cos06 / 2.0) + (
                              c5y * c5z * l5m * vel[0] * sin13 / 2.0) + (cgy * cgz * mg * vel[0] * sin13 / 2.0) - (
                              c5y * l5m * l3z * vel[1] * sin31 / 2.0) - (cgy * l3z * mg * vel[1] * sin31 / 2.0) - (
                              c5x * l5m * l5x * vel[0] * sin13 / 2.0) - (c5y * l5m * l3z * vel[0] * sin27 / 2.0) - (
                              c5y * l5m * l4z * vel[0] * sin25 / 2.0) + (c5y * l5m * l5z * vel[0] * sin13 / 2.0) - (
                              cgx * l5x * mg * vel[0] * sin13 / 2.0) - (cgy * l3z * mg * vel[0] * sin27 / 2.0) - (
                              cgy * l4z * mg * vel[0] * sin25 / 2.0) + (cgy * l5z * mg * vel[0] * sin13 / 2.0) + (
                              c2x * c2y * l2m * vel[1] * cos27) - (c2z * l2m * l2x * vel[0] * cos27) - (
                              c4x * l4m * l3z * vel[0] * cos20) - (c4x * l4m * l4z * vel[0] * cos19) - (
                              c5x * l5m * l4z * vel[1] * cos17 / 2.0) + (c5x * l5m * l4z * vel[2] * cos17 / 2.0) - (
                              c2y * c2z * l2m * vel[1] * sin33) - (cgx * l4z * mg * vel[1] * cos17 / 2.0) + (
                              cgx * l4z * mg * vel[2] * cos17 / 2.0) - (l3m * l2x * l3z * vel[0] * cos27) - (
                              l4m * l2x * l3z * vel[0] * cos27) - (l5m * l2x * l3z * vel[0] * cos27) - (
                              l2x * l3z * mg * vel[0] * cos27) - (l5m * l5x * l3z * vel[0] * cos20) - (
                              l5m * l5x * l4z * vel[0] * cos19) - (c2x * l2m * l2x * vel[0] * sin33) - (
                              c3y * l3m * l3z * vel[1] * sin33) - (c4y * l4m * l3z * vel[1] * sin33) - (
                              l5x * l3z * mg * vel[0] * cos20) - (l5x * l4z * mg * vel[0] * cos19) + (
                              c5x * c5z * l5m * vel[1] * cos01 / 2.0) - (c5x * c5z * l5m * vel[2] * cos01 / 2.0) + (
                              c5x * c5z * l5m * vel[3] * cos01 / 2.0) + (cgx * cgz * mg * vel[1] * cos01 / 2.0) - (
                              cgx * cgz * mg * vel[2] * cos01 / 2.0) + (cgx * cgz * mg * vel[3] * cos01 / 2.0) - (
                              c5y * l5m * l4z * vel[1] * sin22 / 2.0) + (c5y * l5m * l4z * vel[2] * sin22 / 2.0) + (
                              c4z * l4m * l3z * vel[0] * sin19) + (c4z * l4m * l4z * vel[0] * sin20) + (
                              c5z * l5m * l3z * vel[0] * sin19) + (c5z * l5m * l4z * vel[0] * sin20) - (
                              cgy * l4z * mg * vel[1] * sin22 / 2.0) + (cgy * l4z * mg * vel[2] * sin22 / 2.0) + (
                              cgz * l3z * mg * vel[0] * sin19) + (cgz * l4z * mg * vel[0] * sin20)
        C[0, 2] = (I5y * vel[0] * sin10 / 8.0) - (I5x * vel[0] * sin09 / 8.0) - (I5x * vel[0] * sin10 / 8.0) + (
                    I5y * vel[0] * sin09 / 8.0) - (Igx * vel[0] * sin10 / 8.0) - (Igx * vel[0] * sin09 / 8.0) + (
                              Igy * vel[0] * sin10 / 8.0) + (Igy * vel[0] * sin09 / 8.0) - (
                              I4x * vel[0] * sin05 / 2.0) - (I5x * vel[0] * sin05 / 4.0) - (
                              I5y * vel[0] * sin05 / 4.0) + (I4z * vel[0] * sin05 / 2.0) + (
                              I5z * vel[0] * sin05 / 2.0) - (Igx * vel[0] * sin05 / 4.0) - (
                              Igy * vel[0] * sin05 / 4.0) + (Igz * vel[0] * sin05 / 2.0) - (
                              I5x * vel[1] * sin01 / 4.0) + (I5x * vel[1] * sin02 / 4.0) + (
                              I5x * vel[2] * sin01 / 4.0) - (I5x * vel[2] * sin02 / 4.0) - (
                              I5x * vel[3] * sin01 / 4.0) + (I5x * vel[3] * sin02 / 4.0) + (
                              I5x * vel[4] * sin01 / 4.0) + (I5x * vel[4] * sin02 / 4.0) + (
                              I5y * vel[1] * sin01 / 4.0) - (I5y * vel[1] * sin02 / 4.0) - (
                              I5y * vel[2] * sin01 / 4.0) + (I5y * vel[2] * sin02 / 4.0) + (
                              I5y * vel[3] * sin01 / 4.0) - (I5y * vel[3] * sin02 / 4.0) - (
                              I5y * vel[4] * sin01 / 4.0) - (I5y * vel[4] * sin02 / 4.0) - (
                              Igx * vel[1] * sin01 / 4.0) + (Igx * vel[1] * sin02 / 4.0) + (
                              Igx * vel[2] * sin01 / 4.0) - (Igx * vel[2] * sin02 / 4.0) - (
                              Igx * vel[3] * sin01 / 4.0) + (Igx * vel[3] * sin02 / 4.0) + (
                              Igx * vel[4] * sin01 / 4.0) + (Igx * vel[4] * sin02 / 4.0) + (
                              Igy * vel[1] * sin01 / 4.0) - (Igy * vel[1] * sin02 / 4.0) - (
                              Igy * vel[2] * sin01 / 4.0) + (Igy * vel[2] * sin02 / 4.0) + (
                              Igy * vel[3] * sin01 / 4.0) - (Igy * vel[3] * sin02 / 4.0) - (
                              Igy * vel[4] * sin01 / 4.0) - (Igy * vel[4] * sin02 / 4.0) - (
                              I3x * vel[0] * sin18 / 2.0) + (I3z * vel[0] * sin18 / 2.0) + (
                              I5z * vel[4] * sin06 / 2.0) + (Igz * vel[4] * sin06 / 2.0) + (
                              pow(c5x, 2.0) * l5m * vel[0] * sin10 / 8.0) + (
                              pow(c5x, 2.0) * l5m * vel[0] * sin09 / 8.0) - (
                              pow(c5y, 2.0) * l5m * vel[0] * sin10 / 8.0) - (
                              pow(c5y, 2.0) * l5m * vel[0] * sin09 / 8.0) + (
                              pow(cgx, 2.0) * mg * vel[0] * sin10 / 8.0) + (
                              pow(cgx, 2.0) * mg * vel[0] * sin09 / 8.0) - (
                              pow(cgy, 2.0) * mg * vel[0] * sin10 / 8.0) - (
                              pow(cgy, 2.0) * mg * vel[0] * sin09 / 8.0) + (
                              pow(c4x, 2.0) * l4m * vel[0] * sin05 / 2.0) + (
                              pow(c5x, 2.0) * l5m * vel[0] * sin05 / 4.0) + (
                              pow(c5y, 2.0) * l5m * vel[0] * sin05 / 4.0) - (
                              pow(c4z, 2.0) * l4m * vel[0] * sin05 / 2.0) - (
                              pow(c5z, 2.0) * l5m * vel[0] * sin05 / 2.0) + (
                              pow(cgx, 2.0) * mg * vel[0] * sin05 / 4.0) + (
                              pow(cgy, 2.0) * mg * vel[0] * sin05 / 4.0) - (
                              pow(cgz, 2.0) * mg * vel[0] * sin05 / 2.0) + (
                              l5m * pow(l5x, 2.0) * vel[0] * sin05 / 2.0) - (
                              l5m * pow(l5z, 2.0) * vel[0] * sin05 / 2.0) + (
                              pow(l5x, 2.0) * mg * vel[0] * sin05 / 2.0) - (
                              pow(l5z, 2.0) * mg * vel[0] * sin05 / 2.0) + (
                              pow(c5x, 2.0) * l5m * vel[1] * sin01 / 4.0) - (
                              pow(c5x, 2.0) * l5m * vel[1] * sin02 / 4.0) - (
                              pow(c5x, 2.0) * l5m * vel[2] * sin01 / 4.0) + (
                              pow(c5x, 2.0) * l5m * vel[2] * sin02 / 4.0) + (
                              pow(c5x, 2.0) * l5m * vel[3] * sin01 / 4.0) - (
                              pow(c5x, 2.0) * l5m * vel[3] * sin02 / 4.0) - (
                              pow(c5x, 2.0) * l5m * vel[4] * sin01 / 4.0) - (
                              pow(c5x, 2.0) * l5m * vel[4] * sin02 / 4.0) - (
                              pow(c5y, 2.0) * l5m * vel[1] * sin01 / 4.0) + (
                              pow(c5y, 2.0) * l5m * vel[1] * sin02 / 4.0) + (
                              pow(c5y, 2.0) * l5m * vel[2] * sin01 / 4.0) - (
                              pow(c5y, 2.0) * l5m * vel[2] * sin02 / 4.0) - (
                              pow(c5y, 2.0) * l5m * vel[3] * sin01 / 4.0) + (
                              pow(c5y, 2.0) * l5m * vel[3] * sin02 / 4.0) + (
                              pow(c5y, 2.0) * l5m * vel[4] * sin01 / 4.0) + (
                              pow(c5y, 2.0) * l5m * vel[4] * sin02 / 4.0) + (
                              pow(cgx, 2.0) * mg * vel[1] * sin01 / 4.0) - (
                              pow(cgx, 2.0) * mg * vel[1] * sin02 / 4.0) - (
                              pow(cgx, 2.0) * mg * vel[2] * sin01 / 4.0) + (
                              pow(cgx, 2.0) * mg * vel[2] * sin02 / 4.0) + (
                              pow(cgx, 2.0) * mg * vel[3] * sin01 / 4.0) - (
                              pow(cgx, 2.0) * mg * vel[3] * sin02 / 4.0) - (
                              pow(cgx, 2.0) * mg * vel[4] * sin01 / 4.0) - (
                              pow(cgx, 2.0) * mg * vel[4] * sin02 / 4.0) - (
                              pow(cgy, 2.0) * mg * vel[1] * sin01 / 4.0) + (
                              pow(cgy, 2.0) * mg * vel[1] * sin02 / 4.0) + (
                              pow(cgy, 2.0) * mg * vel[2] * sin01 / 4.0) - (
                              pow(cgy, 2.0) * mg * vel[2] * sin02 / 4.0) - (
                              pow(cgy, 2.0) * mg * vel[3] * sin01 / 4.0) + (
                              pow(cgy, 2.0) * mg * vel[3] * sin02 / 4.0) + (
                              pow(cgy, 2.0) * mg * vel[4] * sin01 / 4.0) + (
                              pow(cgy, 2.0) * mg * vel[4] * sin02 / 4.0) + (
                              pow(c3x, 2.0) * l3m * vel[0] * sin18 / 2.0) - (
                              pow(c3z, 2.0) * l3m * vel[0] * sin18 / 2.0) - (
                              l4m * pow(l4z, 2.0) * vel[0] * sin18 / 2.0) - (
                              l5m * pow(l4z, 2.0) * vel[0] * sin18 / 2.0) - (
                              pow(l4z, 2.0) * mg * vel[0] * sin18 / 2.0) + (
                              pow(c5x, 2.0) * l5m * vel[4] * sin06 / 2.0) + (
                              pow(c5y, 2.0) * l5m * vel[4] * sin06 / 2.0) + (
                              pow(cgx, 2.0) * mg * vel[4] * sin06 / 2.0) + (
                              pow(cgy, 2.0) * mg * vel[4] * sin06 / 2.0) - (l3z * l4z * mg * vel[0] * sin24 / 2.0) + (
                              c5y * l5m * l2x * vel[0] * cos01 / 2.0) - (c5y * l5m * l5x * vel[1] * cos01 / 2.0) + (
                              c5y * l5m * l5x * vel[2] * cos01 / 2.0) - (c5y * l5m * l5x * vel[3] * cos01 / 2.0) - (
                              c5x * l5m * l5z * vel[1] * cos01 / 2.0) + (c5x * l5m * l5z * vel[2] * cos01 / 2.0) - (
                              c5x * l5m * l5z * vel[3] * cos01 / 2.0) + (cgy * l2x * mg * vel[0] * cos01 / 2.0) - (
                              cgy * l5x * mg * vel[1] * cos01 / 2.0) + (cgy * l5x * mg * vel[2] * cos01 / 2.0) - (
                              cgy * l5x * mg * vel[3] * cos01 / 2.0) - (cgx * l5z * mg * vel[1] * cos01 / 2.0) + (
                              cgx * l5z * mg * vel[2] * cos01 / 2.0) - (cgx * l5z * mg * vel[3] * cos01 / 2.0) - (
                              l5m * l3z * l5z * vel[0] * sin19 / 2.0) - (l5m * l4z * l5z * vel[0] * sin20) - (
                              l3z * l5z * mg * vel[0] * sin19 / 2.0) - (l4z * l5z * mg * vel[0] * sin20) + (
                              c5y * c5z * l5m * vel[1] * sin03 / 2.0) - (c5y * c5z * l5m * vel[2] * sin03 / 2.0) + (
                              c5y * c5z * l5m * vel[3] * sin03 / 2.0) + (cgy * cgz * mg * vel[1] * sin03 / 2.0) - (
                              cgy * cgz * mg * vel[2] * sin03 / 2.0) + (cgy * cgz * mg * vel[3] * sin03 / 2.0) - (
                              c3x * c3y * l3m * vel[1] * cos18) + (c3x * c3y * l3m * vel[2] * cos18) + (
                              c5x * l5m * l2x * vel[0] * sin03 / 2.0) - (c5x * l5m * l5x * vel[1] * sin03 / 2.0) + (
                              c5x * l5m * l5x * vel[2] * sin03 / 2.0) - (c5x * l5m * l5x * vel[3] * sin03 / 2.0) + (
                              c5y * l5m * l5z * vel[1] * sin03 / 2.0) - (c5y * l5m * l5z * vel[2] * sin03 / 2.0) + (
                              c5y * l5m * l5z * vel[3] * sin03 / 2.0) + (cgx * l2x * mg * vel[0] * sin03 / 2.0) - (
                              cgx * l5x * mg * vel[1] * sin03 / 2.0) + (cgx * l5x * mg * vel[2] * sin03 / 2.0) - (
                              cgx * l5x * mg * vel[3] * sin03 / 2.0) + (cgy * l5z * mg * vel[1] * sin03 / 2.0) - (
                              cgy * l5z * mg * vel[2] * sin03 / 2.0) + (cgy * l5z * mg * vel[3] * sin03 / 2.0) - (
                              c5x * c5y * l5m * vel[0] * cos23 / 4.0) + (c5x * c5y * l5m * vel[0] * cos26 / 4.0) + (
                              c5x * c5z * l5m * vel[0] * cos05 / 2.0) - (cgx * cgy * mg * vel[0] * cos23 / 4.0) + (
                              cgx * cgy * mg * vel[0] * cos26 / 4.0) + (cgx * cgz * mg * vel[0] * cos05 / 2.0) + (
                              c3z * l3m * l2x * vel[0] * cos18) - (c4x * l4m * l3z * vel[0] * cos12 / 2.0) - (
                              c5y * l5m * l5x * vel[0] * cos05 / 2.0) + (c5x * l5m * l5z * vel[0] * cos05 / 2.0) + (
                              c3y * c3z * l3m * vel[1] * sin30) - (c3y * c3z * l3m * vel[2] * sin30) - (
                              cgy * l5x * mg * vel[0] * cos05 / 2.0) + (cgx * l5z * mg * vel[0] * cos05 / 2.0) + (
                              l4m * l2x * l4z * vel[0] * cos18) + (l5m * l2x * l4z * vel[0] * cos18) - (
                              l5m * l5x * l3z * vel[0] * cos12 / 2.0) + (l2x * l4z * mg * vel[0] * cos18) - (
                              l5x * l3z * mg * vel[0] * cos12 / 2.0) + (c5y * c5z * l5m * vel[0] * sin12 / 2.0) + (
                              cgy * cgz * mg * vel[0] * sin12 / 2.0) + (c3x * l3m * l2x * vel[0] * sin30) + (
                              c4y * l4m * l4z * vel[1] * sin30) - (c4y * l4m * l4z * vel[2] * sin30) - (
                              c4z * l4m * l3z * vel[0] * sin11 / 2.0) - (c5z * l5m * l3z * vel[0] * sin11 / 2.0) - (
                              cgz * l3z * mg * vel[0] * sin11 / 2.0) + (c5x * l5m * l5x * vel[0] * sin12 / 2.0) + (
                              c5y * l5m * l5z * vel[0] * sin12 / 2.0) + (cgx * l5x * mg * vel[0] * sin12 / 2.0) + (
                              cgy * l5z * mg * vel[0] * sin12 / 2.0) - (l5m * l3z * l5z * vel[0] * sin11 / 2.0) - (
                              l3z * l5z * mg * vel[0] * sin11 / 2.0) + (c4x * c4z * l4m * vel[0] * cos15) + (
                              c5z * l5m * l5x * vel[0] * cos15) + (cgz * l5x * mg * vel[0] * cos15) + (
                              l5m * l5x * l5z * vel[0] * cos15) + (l5x * l5z * mg * vel[0] * cos15) - (
                              c5x * c5y * l5m * vel[1] * cos08 / 2.0) - (c5x * c5y * l5m * vel[1] * cos10 / 2.0) + (
                              c5x * c5y * l5m * vel[2] * cos08 / 2.0) + (c5x * c5y * l5m * vel[2] * cos10 / 2.0) - (
                              c5x * c5y * l5m * vel[3] * cos08 / 2.0) - (c5x * c5y * l5m * vel[3] * cos10 / 2.0) + (
                              c5x * c5y * l5m * vel[4] * cos08 / 2.0) - (c5x * c5y * l5m * vel[4] * cos10 / 2.0) + (
                              c5x * c5z * l5m * vel[1] * cos00 / 2.0) - (c5x * c5z * l5m * vel[2] * cos00 / 2.0) + (
                              c5x * c5z * l5m * vel[3] * cos00 / 2.0) - (cgx * cgy * mg * vel[1] * cos08 / 2.0) - (
                              cgx * cgy * mg * vel[1] * cos10 / 2.0) + (cgx * cgy * mg * vel[2] * cos08 / 2.0) + (
                              cgx * cgy * mg * vel[2] * cos10 / 2.0) - (cgx * cgy * mg * vel[3] * cos08 / 2.0) - (
                              cgx * cgy * mg * vel[3] * cos10 / 2.0) + (cgx * cgy * mg * vel[4] * cos08 / 2.0) - (
                              cgx * cgy * mg * vel[4] * cos10 / 2.0) + (cgx * cgz * mg * vel[1] * cos00 / 2.0) - (
                              cgx * cgz * mg * vel[2] * cos00 / 2.0) + (cgx * cgz * mg * vel[3] * cos00 / 2.0) - (
                              c5z * l5m * l5z * vel[0] * sin05) - (cgz * l5z * mg * vel[0] * sin05) - (
                              c5y * l5m * l2x * vel[0] * cos00 / 2.0) - (c5y * l5m * l5x * vel[1] * cos00 / 2.0) + (
                              c5y * l5m * l5x * vel[2] * cos00 / 2.0) - (c5y * l5m * l5x * vel[3] * cos00 / 2.0) + (
                              c5x * l5m * l3z * vel[0] * cos25 / 4.0) + (c5x * l5m * l4z * vel[0] * cos24 / 2.0) + (
                              c5x * l5m * l5z * vel[1] * cos00 / 2.0) - (c5x * l5m * l5z * vel[2] * cos00 / 2.0) + (
                              c5x * l5m * l5z * vel[3] * cos00 / 2.0) - (cgy * l2x * mg * vel[0] * cos00 / 2.0) - (
                              cgy * l5x * mg * vel[1] * cos00 / 2.0) + (cgy * l5x * mg * vel[2] * cos00 / 2.0) - (
                              cgy * l5x * mg * vel[3] * cos00 / 2.0) + (cgx * l3z * mg * vel[0] * cos25 / 4.0) + (
                              cgx * l4z * mg * vel[0] * cos24 / 2.0) + (cgx * l5z * mg * vel[1] * cos00 / 2.0) - (
                              cgx * l5z * mg * vel[2] * cos00 / 2.0) + (cgx * l5z * mg * vel[3] * cos00 / 2.0) + (
                              c5y * c5z * l5m * vel[1] * sin04 / 2.0) - (c5y * c5z * l5m * vel[2] * sin04 / 2.0) + (
                              c5y * c5z * l5m * vel[3] * sin04 / 2.0) + (cgy * cgz * mg * vel[1] * sin04 / 2.0) - (
                              cgy * cgz * mg * vel[2] * sin04 / 2.0) + (cgy * cgz * mg * vel[3] * sin04 / 2.0) + (
                              c3x * c3z * l3m * vel[0] * cos31) + (c5x * l5m * l2x * vel[0] * sin04 / 2.0) + (
                              c5x * l5m * l5x * vel[1] * sin04 / 2.0) - (c5x * l5m * l5x * vel[2] * sin04 / 2.0) + (
                              c5x * l5m * l5x * vel[3] * sin04 / 2.0) - (c5y * l5m * l3z * vel[0] * sin26 / 4.0) - (
                              c5y * l5m * l4z * vel[0] * sin28 / 2.0) + (c5y * l5m * l5z * vel[1] * sin04 / 2.0) - (
                              c5y * l5m * l5z * vel[2] * sin04 / 2.0) + (c5y * l5m * l5z * vel[3] * sin04 / 2.0) + (
                              cgx * l2x * mg * vel[0] * sin04 / 2.0) + (cgx * l5x * mg * vel[1] * sin04 / 2.0) - (
                              cgx * l5x * mg * vel[2] * sin04 / 2.0) + (cgx * l5x * mg * vel[3] * sin04 / 2.0) - (
                              cgy * l3z * mg * vel[0] * sin26 / 4.0) - (cgy * l4z * mg * vel[0] * sin28 / 2.0) + (
                              cgy * l5z * mg * vel[1] * sin04 / 2.0) - (cgy * l5z * mg * vel[2] * sin04 / 2.0) + (
                              cgy * l5z * mg * vel[3] * sin04 / 2.0) + (c3x * l3m * l3z * vel[0] * cos32 / 2.0) - (
                              c4x * c4y * l4m * vel[1] * cos04) + (c4x * c4y * l4m * vel[2] * cos04) - (
                              c4x * c4y * l4m * vel[3] * cos04) - (c3z * l3m * l3z * vel[0] * sin29 / 2.0) + (
                              c4z * l4m * l2x * vel[0] * cos04) + (c5z * l5m * l2x * vel[0] * cos04) - (
                              c5x * l5m * l3z * vel[0] * cos14 / 4.0) - (c5x * l5m * l4z * vel[1] * cos16 / 2.0) + (
                              c5x * l5m * l4z * vel[2] * cos16 / 2.0) + (cgz * l2x * mg * vel[0] * cos04) - (
                              cgx * l3z * mg * vel[0] * cos14 / 4.0) - (cgx * l4z * mg * vel[1] * cos16 / 2.0) + (
                              cgx * l4z * mg * vel[2] * cos16 / 2.0) - (l4m * l3z * l4z * vel[0] * sin29 / 2.0) - (
                              l5m * l3z * l4z * vel[0] * sin29 / 2.0) - (l3z * l4z * mg * vel[0] * sin29 / 2.0) + (
                              c4y * c4z * l4m * vel[1] * sin06) - (c4y * c4z * l4m * vel[2] * sin06) + (
                              c4y * c4z * l4m * vel[3] * sin06) + (l5m * l2x * l5z * vel[0] * cos04) + (
                              l2x * l5z * mg * vel[0] * cos04) + (c4x * l4m * l2x * vel[0] * sin06) + (
                              c5y * l5m * l3z * vel[0] * sin16 / 4.0) + (c5y * l5m * l4z * vel[1] * sin21 / 2.0) - (
                              c5y * l5m * l4z * vel[2] * sin21 / 2.0) + (cgy * l3z * mg * vel[0] * sin16 / 4.0) + (
                              cgy * l4z * mg * vel[1] * sin21 / 2.0) - (cgy * l4z * mg * vel[2] * sin21 / 2.0) + (
                              l5m * l2x * l5x * vel[0] * sin06) + (l2x * l5x * mg * vel[0] * sin06) + (
                              c5x * c5z * l5m * vel[0] * cos06 / 2.0) + (cgx * cgz * mg * vel[0] * cos06 / 2.0) + (
                              c5y * l5m * l5x * vel[0] * cos06 / 2.0) + (c5x * l5m * l3z * vel[0] * cos22 / 4.0) + (
                              c5x * l5m * l4z * vel[0] * cos21 / 2.0) + (c5x * l5m * l5z * vel[0] * cos06 / 2.0) + (
                              cgy * l5x * mg * vel[0] * cos06 / 2.0) + (cgx * l3z * mg * vel[0] * cos22 / 4.0) + (
                              cgx * l4z * mg * vel[0] * cos21 / 2.0) + (cgx * l5z * mg * vel[0] * cos06 / 2.0) - (
                              c5y * c5z * l5m * vel[0] * sin13 / 2.0) - (cgy * cgz * mg * vel[0] * sin13 / 2.0) + (
                              c5x * l5m * l5x * vel[0] * sin13 / 2.0) + (c5y * l5m * l3z * vel[0] * sin27 / 4.0) + (
                              c5y * l5m * l4z * vel[0] * sin25 / 2.0) - (c5y * l5m * l5z * vel[0] * sin13 / 2.0) + (
                              cgx * l5x * mg * vel[0] * sin13 / 2.0) + (cgy * l3z * mg * vel[0] * sin27 / 4.0) + (
                              cgy * l4z * mg * vel[0] * sin25 / 2.0) - (cgy * l5z * mg * vel[0] * sin13 / 2.0) - (
                              c3x * l3m * l3z * vel[0] * cos30 / 2.0) + (c4x * l4m * l3z * vel[0] * cos20 / 2.0) + (
                              c4x * l4m * l4z * vel[0] * cos19) - (c5x * l5m * l3z * vel[0] * cos13 / 4.0) + (
                              c5x * l5m * l4z * vel[1] * cos17 / 2.0) - (c5x * l5m * l4z * vel[2] * cos17 / 2.0) - (
                              cgx * l3z * mg * vel[0] * cos13 / 4.0) + (cgx * l4z * mg * vel[1] * cos17 / 2.0) - (
                              cgx * l4z * mg * vel[2] * cos17 / 2.0) + (l5m * l5x * l3z * vel[0] * cos20 / 2.0) + (
                              l5m * l5x * l4z * vel[0] * cos19) - (c3z * l3m * l3z * vel[0] * sin24 / 2.0) + (
                              l5x * l3z * mg * vel[0] * cos20 / 2.0) + (l5x * l4z * mg * vel[0] * cos19) - (
                              c5x * c5z * l5m * vel[1] * cos01 / 2.0) + (c5x * c5z * l5m * vel[2] * cos01 / 2.0) - (
                              c5x * c5z * l5m * vel[3] * cos01 / 2.0) - (cgx * cgz * mg * vel[1] * cos01 / 2.0) + (
                              cgx * cgz * mg * vel[2] * cos01 / 2.0) - (cgx * cgz * mg * vel[3] * cos01 / 2.0) - (
                              c5y * l5m * l3z * vel[0] * sin17 / 4.0) + (c5y * l5m * l4z * vel[1] * sin22 / 2.0) - (
                              c5y * l5m * l4z * vel[2] * sin22 / 2.0) - (c4z * l4m * l3z * vel[0] * sin19 / 2.0) - (
                              c4z * l4m * l4z * vel[0] * sin20) - (c5z * l5m * l3z * vel[0] * sin19 / 2.0) - (
                              c5z * l5m * l4z * vel[0] * sin20) - (cgy * l3z * mg * vel[0] * sin17 / 4.0) + (
                              cgy * l4z * mg * vel[1] * sin22 / 2.0) - (cgy * l4z * mg * vel[2] * sin22 / 2.0) - (
                              cgz * l3z * mg * vel[0] * sin19 / 2.0) - (cgz * l4z * mg * vel[0] * sin20) - (
                              l4m * l3z * l4z * vel[0] * sin24 / 2.0) - (l5m * l3z * l4z * vel[0] * sin24 / 2.0)
        C[0, 3] = (I5x * vel[0] * sin10 / 8.0) + (I5x * vel[0] * sin09 / 8.0) - (I5y * vel[0] * sin10 / 8.0) - (
                    I5y * vel[0] * sin09 / 8.0) + (Igx * vel[0] * sin10 / 8.0) + (Igx * vel[0] * sin09 / 8.0) - (
                              Igy * vel[0] * sin10 / 8.0) - (Igy * vel[0] * sin09 / 8.0) + (
                              I4x * vel[0] * sin05 / 2.0) + (I5x * vel[0] * sin05 / 4.0) + (
                              I5y * vel[0] * sin05 / 4.0) - (I4z * vel[0] * sin05 / 2.0) - (
                              I5z * vel[0] * sin05 / 2.0) + (Igx * vel[0] * sin05 / 4.0) + (
                              Igy * vel[0] * sin05 / 4.0) - (Igz * vel[0] * sin05 / 2.0) + (
                              I5x * vel[1] * sin01 / 4.0) - (I5x * vel[1] * sin02 / 4.0) - (
                              I5x * vel[2] * sin01 / 4.0) + (I5x * vel[2] * sin02 / 4.0) + (
                              I5x * vel[3] * sin01 / 4.0) - (I5x * vel[3] * sin02 / 4.0) - (
                              I5x * vel[4] * sin01 / 4.0) - (I5x * vel[4] * sin02 / 4.0) - (
                              I5y * vel[1] * sin01 / 4.0) + (I5y * vel[1] * sin02 / 4.0) + (
                              I5y * vel[2] * sin01 / 4.0) - (I5y * vel[2] * sin02 / 4.0) - (
                              I5y * vel[3] * sin01 / 4.0) + (I5y * vel[3] * sin02 / 4.0) + (
                              I5y * vel[4] * sin01 / 4.0) + (I5y * vel[4] * sin02 / 4.0) + (
                              Igx * vel[1] * sin01 / 4.0) - (Igx * vel[1] * sin02 / 4.0) - (
                              Igx * vel[2] * sin01 / 4.0) + (Igx * vel[2] * sin02 / 4.0) + (
                              Igx * vel[3] * sin01 / 4.0) - (Igx * vel[3] * sin02 / 4.0) - (
                              Igx * vel[4] * sin01 / 4.0) - (Igx * vel[4] * sin02 / 4.0) - (
                              Igy * vel[1] * sin01 / 4.0) + (Igy * vel[1] * sin02 / 4.0) + (
                              Igy * vel[2] * sin01 / 4.0) - (Igy * vel[2] * sin02 / 4.0) - (
                              Igy * vel[3] * sin01 / 4.0) + (Igy * vel[3] * sin02 / 4.0) + (
                              Igy * vel[4] * sin01 / 4.0) + (Igy * vel[4] * sin02 / 4.0) - (
                              I5z * vel[4] * sin06 / 2.0) - (Igz * vel[4] * sin06 / 2.0) - (
                              pow(c5x, 2.0) * l5m * vel[0] * sin10 / 8.0) - (
                              pow(c5x, 2.0) * l5m * vel[0] * sin09 / 8.0) + (
                              pow(c5y, 2.0) * l5m * vel[0] * sin10 / 8.0) + (
                              pow(c5y, 2.0) * l5m * vel[0] * sin09 / 8.0) - (
                              pow(cgx, 2.0) * mg * vel[0] * sin10 / 8.0) - (
                              pow(cgx, 2.0) * mg * vel[0] * sin09 / 8.0) + (
                              pow(cgy, 2.0) * mg * vel[0] * sin10 / 8.0) + (
                              pow(cgy, 2.0) * mg * vel[0] * sin09 / 8.0) - (
                              pow(c4x, 2.0) * l4m * vel[0] * sin05 / 2.0) - (
                              pow(c5x, 2.0) * l5m * vel[0] * sin05 / 4.0) - (
                              pow(c5y, 2.0) * l5m * vel[0] * sin05 / 4.0) + (
                              pow(c4z, 2.0) * l4m * vel[0] * sin05 / 2.0) + (
                              pow(c5z, 2.0) * l5m * vel[0] * sin05 / 2.0) - (
                              pow(cgx, 2.0) * mg * vel[0] * sin05 / 4.0) - (
                              pow(cgy, 2.0) * mg * vel[0] * sin05 / 4.0) + (
                              pow(cgz, 2.0) * mg * vel[0] * sin05 / 2.0) - (
                              l5m * pow(l5x, 2.0) * vel[0] * sin05 / 2.0) + (
                              l5m * pow(l5z, 2.0) * vel[0] * sin05 / 2.0) - (
                              pow(l5x, 2.0) * mg * vel[0] * sin05 / 2.0) + (
                              pow(l5z, 2.0) * mg * vel[0] * sin05 / 2.0) - (
                              pow(c5x, 2.0) * l5m * vel[1] * sin01 / 4.0) + (
                              pow(c5x, 2.0) * l5m * vel[1] * sin02 / 4.0) + (
                              pow(c5x, 2.0) * l5m * vel[2] * sin01 / 4.0) - (
                              pow(c5x, 2.0) * l5m * vel[2] * sin02 / 4.0) - (
                              pow(c5x, 2.0) * l5m * vel[3] * sin01 / 4.0) + (
                              pow(c5x, 2.0) * l5m * vel[3] * sin02 / 4.0) + (
                              pow(c5x, 2.0) * l5m * vel[4] * sin01 / 4.0) + (
                              pow(c5x, 2.0) * l5m * vel[4] * sin02 / 4.0) + (
                              pow(c5y, 2.0) * l5m * vel[1] * sin01 / 4.0) - (
                              pow(c5y, 2.0) * l5m * vel[1] * sin02 / 4.0) - (
                              pow(c5y, 2.0) * l5m * vel[2] * sin01 / 4.0) + (
                              pow(c5y, 2.0) * l5m * vel[2] * sin02 / 4.0) + (
                              pow(c5y, 2.0) * l5m * vel[3] * sin01 / 4.0) - (
                              pow(c5y, 2.0) * l5m * vel[3] * sin02 / 4.0) - (
                              pow(c5y, 2.0) * l5m * vel[4] * sin01 / 4.0) - (
                              pow(c5y, 2.0) * l5m * vel[4] * sin02 / 4.0) - (
                              pow(cgx, 2.0) * mg * vel[1] * sin01 / 4.0) + (
                              pow(cgx, 2.0) * mg * vel[1] * sin02 / 4.0) + (
                              pow(cgx, 2.0) * mg * vel[2] * sin01 / 4.0) - (
                              pow(cgx, 2.0) * mg * vel[2] * sin02 / 4.0) - (
                              pow(cgx, 2.0) * mg * vel[3] * sin01 / 4.0) + (
                              pow(cgx, 2.0) * mg * vel[3] * sin02 / 4.0) + (
                              pow(cgx, 2.0) * mg * vel[4] * sin01 / 4.0) + (
                              pow(cgx, 2.0) * mg * vel[4] * sin02 / 4.0) + (
                              pow(cgy, 2.0) * mg * vel[1] * sin01 / 4.0) - (
                              pow(cgy, 2.0) * mg * vel[1] * sin02 / 4.0) - (
                              pow(cgy, 2.0) * mg * vel[2] * sin01 / 4.0) + (
                              pow(cgy, 2.0) * mg * vel[2] * sin02 / 4.0) + (
                              pow(cgy, 2.0) * mg * vel[3] * sin01 / 4.0) - (
                              pow(cgy, 2.0) * mg * vel[3] * sin02 / 4.0) - (
                              pow(cgy, 2.0) * mg * vel[4] * sin01 / 4.0) - (
                              pow(cgy, 2.0) * mg * vel[4] * sin02 / 4.0) - (
                              pow(c5x, 2.0) * l5m * vel[4] * sin06 / 2.0) - (
                              pow(c5y, 2.0) * l5m * vel[4] * sin06 / 2.0) - (
                              pow(cgx, 2.0) * mg * vel[4] * sin06 / 2.0) - (
                              pow(cgy, 2.0) * mg * vel[4] * sin06 / 2.0) - (l4z * l5z * mg * vel[0] * sin08 / 2.0) - (
                              c5y * l5m * l2x * vel[0] * cos01 / 2.0) + (c5y * l5m * l5x * vel[1] * cos01 / 2.0) - (
                              c5y * l5m * l5x * vel[2] * cos01 / 2.0) + (c5y * l5m * l5x * vel[3] * cos01 / 2.0) + (
                              c5x * l5m * l5z * vel[1] * cos01 / 2.0) - (c5x * l5m * l5z * vel[2] * cos01 / 2.0) + (
                              c5x * l5m * l5z * vel[3] * cos01 / 2.0) - (cgy * l2x * mg * vel[0] * cos01 / 2.0) + (
                              cgy * l5x * mg * vel[1] * cos01 / 2.0) - (cgy * l5x * mg * vel[2] * cos01 / 2.0) + (
                              cgy * l5x * mg * vel[3] * cos01 / 2.0) + (cgx * l5z * mg * vel[1] * cos01 / 2.0) - (
                              cgx * l5z * mg * vel[2] * cos01 / 2.0) + (cgx * l5z * mg * vel[3] * cos01 / 2.0) + (
                              l5m * l3z * l5z * vel[0] * sin19 / 2.0) + (l5m * l4z * l5z * vel[0] * sin20 / 2.0) + (
                              l3z * l5z * mg * vel[0] * sin19 / 2.0) + (l4z * l5z * mg * vel[0] * sin20 / 2.0) - (
                              c5y * c5z * l5m * vel[1] * sin03 / 2.0) + (c5y * c5z * l5m * vel[2] * sin03 / 2.0) - (
                              c5y * c5z * l5m * vel[3] * sin03 / 2.0) - (cgy * cgz * mg * vel[1] * sin03 / 2.0) + (
                              cgy * cgz * mg * vel[2] * sin03 / 2.0) - (cgy * cgz * mg * vel[3] * sin03 / 2.0) - (
                              c5x * l5m * l2x * vel[0] * sin03 / 2.0) + (c5x * l5m * l5x * vel[1] * sin03 / 2.0) - (
                              c5x * l5m * l5x * vel[2] * sin03 / 2.0) + (c5x * l5m * l5x * vel[3] * sin03 / 2.0) - (
                              c5y * l5m * l5z * vel[1] * sin03 / 2.0) + (c5y * l5m * l5z * vel[2] * sin03 / 2.0) - (
                              c5y * l5m * l5z * vel[3] * sin03 / 2.0) - (cgx * l2x * mg * vel[0] * sin03 / 2.0) + (
                              cgx * l5x * mg * vel[1] * sin03 / 2.0) - (cgx * l5x * mg * vel[2] * sin03 / 2.0) + (
                              cgx * l5x * mg * vel[3] * sin03 / 2.0) - (cgy * l5z * mg * vel[1] * sin03 / 2.0) + (
                              cgy * l5z * mg * vel[2] * sin03 / 2.0) - (cgy * l5z * mg * vel[3] * sin03 / 2.0) + (
                              c5x * c5y * l5m * vel[0] * cos23 / 4.0) - (c5x * c5y * l5m * vel[0] * cos26 / 4.0) - (
                              c5x * c5z * l5m * vel[0] * cos05 / 2.0) + (cgx * cgy * mg * vel[0] * cos23 / 4.0) - (
                              cgx * cgy * mg * vel[0] * cos26 / 4.0) - (cgx * cgz * mg * vel[0] * cos05 / 2.0) + (
                              c4x * l4m * l3z * vel[0] * cos12 / 2.0) + (c5x * l5m * l4z * vel[0] * cos09 / 4.0) + (
                              cgx * l4z * mg * vel[0] * cos09 / 4.0) + (c5y * l5m * l5x * vel[0] * cos05 / 2.0) - (
                              c5x * l5m * l5z * vel[0] * cos05 / 2.0) + (cgy * l5x * mg * vel[0] * cos05 / 2.0) - (
                              cgx * l5z * mg * vel[0] * cos05 / 2.0) + (l5m * l5x * l3z * vel[0] * cos12 / 2.0) + (
                              l5x * l3z * mg * vel[0] * cos12 / 2.0) - (c5y * c5z * l5m * vel[0] * sin12 / 2.0) - (
                              cgy * cgz * mg * vel[0] * sin12 / 2.0) + (c5y * l5m * l4z * vel[0] * sin15 / 4.0) + (
                              c4z * l4m * l3z * vel[0] * sin11 / 2.0) + (c5z * l5m * l3z * vel[0] * sin11 / 2.0) + (
                              cgy * l4z * mg * vel[0] * sin15 / 4.0) + (cgz * l3z * mg * vel[0] * sin11 / 2.0) - (
                              c5x * l5m * l5x * vel[0] * sin12 / 2.0) - (c5y * l5m * l5z * vel[0] * sin12 / 2.0) - (
                              cgx * l5x * mg * vel[0] * sin12 / 2.0) - (cgy * l5z * mg * vel[0] * sin12 / 2.0) + (
                              l5m * l3z * l5z * vel[0] * sin11 / 2.0) + (l3z * l5z * mg * vel[0] * sin11 / 2.0) - (
                              c4x * c4z * l4m * vel[0] * cos15) - (c5z * l5m * l5x * vel[0] * cos15) - (
                              cgz * l5x * mg * vel[0] * cos15) - (l5m * l5x * l5z * vel[0] * cos15) - (
                              l5x * l5z * mg * vel[0] * cos15) + (c5x * c5y * l5m * vel[1] * cos08 / 2.0) + (
                              c5x * c5y * l5m * vel[1] * cos10 / 2.0) - (c5x * c5y * l5m * vel[2] * cos08 / 2.0) - (
                              c5x * c5y * l5m * vel[2] * cos10 / 2.0) + (c5x * c5y * l5m * vel[3] * cos08 / 2.0) + (
                              c5x * c5y * l5m * vel[3] * cos10 / 2.0) - (c5x * c5y * l5m * vel[4] * cos08 / 2.0) + (
                              c5x * c5y * l5m * vel[4] * cos10 / 2.0) - (c5x * c5z * l5m * vel[1] * cos00 / 2.0) + (
                              c5x * c5z * l5m * vel[2] * cos00 / 2.0) - (c5x * c5z * l5m * vel[3] * cos00 / 2.0) + (
                              cgx * cgy * mg * vel[1] * cos08 / 2.0) + (cgx * cgy * mg * vel[1] * cos10 / 2.0) - (
                              cgx * cgy * mg * vel[2] * cos08 / 2.0) - (cgx * cgy * mg * vel[2] * cos10 / 2.0) + (
                              cgx * cgy * mg * vel[3] * cos08 / 2.0) + (cgx * cgy * mg * vel[3] * cos10 / 2.0) - (
                              cgx * cgy * mg * vel[4] * cos08 / 2.0) + (cgx * cgy * mg * vel[4] * cos10 / 2.0) - (
                              cgx * cgz * mg * vel[1] * cos00 / 2.0) + (cgx * cgz * mg * vel[2] * cos00 / 2.0) - (
                              cgx * cgz * mg * vel[3] * cos00 / 2.0) + (c5z * l5m * l5z * vel[0] * sin05) + (
                              cgz * l5z * mg * vel[0] * sin05) + (c5y * l5m * l2x * vel[0] * cos00 / 2.0) + (
                              c5y * l5m * l5x * vel[1] * cos00 / 2.0) - (c5y * l5m * l5x * vel[2] * cos00 / 2.0) + (
                              c5y * l5m * l5x * vel[3] * cos00 / 2.0) - (c5x * l5m * l3z * vel[0] * cos25 / 4.0) - (
                              c5x * l5m * l4z * vel[0] * cos24 / 4.0) - (c5x * l5m * l5z * vel[1] * cos00 / 2.0) + (
                              c5x * l5m * l5z * vel[2] * cos00 / 2.0) - (c5x * l5m * l5z * vel[3] * cos00 / 2.0) + (
                              cgy * l2x * mg * vel[0] * cos00 / 2.0) + (cgy * l5x * mg * vel[1] * cos00 / 2.0) - (
                              cgy * l5x * mg * vel[2] * cos00 / 2.0) + (cgy * l5x * mg * vel[3] * cos00 / 2.0) - (
                              cgx * l3z * mg * vel[0] * cos25 / 4.0) - (cgx * l4z * mg * vel[0] * cos24 / 4.0) - (
                              cgx * l5z * mg * vel[1] * cos00 / 2.0) + (cgx * l5z * mg * vel[2] * cos00 / 2.0) - (
                              cgx * l5z * mg * vel[3] * cos00 / 2.0) - (c5y * c5z * l5m * vel[1] * sin04 / 2.0) + (
                              c5y * c5z * l5m * vel[2] * sin04 / 2.0) - (c5y * c5z * l5m * vel[3] * sin04 / 2.0) - (
                              cgy * cgz * mg * vel[1] * sin04 / 2.0) + (cgy * cgz * mg * vel[2] * sin04 / 2.0) - (
                              cgy * cgz * mg * vel[3] * sin04 / 2.0) - (c5x * l5m * l2x * vel[0] * sin04 / 2.0) - (
                              c5x * l5m * l5x * vel[1] * sin04 / 2.0) + (c5x * l5m * l5x * vel[2] * sin04 / 2.0) - (
                              c5x * l5m * l5x * vel[3] * sin04 / 2.0) + (c5y * l5m * l3z * vel[0] * sin26 / 4.0) + (
                              c5y * l5m * l4z * vel[0] * sin28 / 4.0) - (c5y * l5m * l5z * vel[1] * sin04 / 2.0) + (
                              c5y * l5m * l5z * vel[2] * sin04 / 2.0) - (c5y * l5m * l5z * vel[3] * sin04 / 2.0) - (
                              cgx * l2x * mg * vel[0] * sin04 / 2.0) - (cgx * l5x * mg * vel[1] * sin04 / 2.0) + (
                              cgx * l5x * mg * vel[2] * sin04 / 2.0) - (cgx * l5x * mg * vel[3] * sin04 / 2.0) + (
                              cgy * l3z * mg * vel[0] * sin26 / 4.0) + (cgy * l4z * mg * vel[0] * sin28 / 4.0) - (
                              cgy * l5z * mg * vel[1] * sin04 / 2.0) + (cgy * l5z * mg * vel[2] * sin04 / 2.0) - (
                              cgy * l5z * mg * vel[3] * sin04 / 2.0) + (c4x * c4y * l4m * vel[1] * cos04) - (
                              c4x * c4y * l4m * vel[2] * cos04) + (c4x * c4y * l4m * vel[3] * cos04) - (
                              c4z * l4m * l2x * vel[0] * cos04) - (c5z * l5m * l2x * vel[0] * cos04) + (
                              c5x * l5m * l3z * vel[0] * cos14 / 4.0) - (cgz * l2x * mg * vel[0] * cos04) + (
                              cgx * l3z * mg * vel[0] * cos14 / 4.0) - (c4y * c4z * l4m * vel[1] * sin06) + (
                              c4y * c4z * l4m * vel[2] * sin06) - (c4y * c4z * l4m * vel[3] * sin06) - (
                              l5m * l2x * l5z * vel[0] * cos04) - (l2x * l5z * mg * vel[0] * cos04) - (
                              c4x * l4m * l2x * vel[0] * sin06) - (c5y * l5m * l3z * vel[0] * sin16 / 4.0) - (
                              cgy * l3z * mg * vel[0] * sin16 / 4.0) - (l5m * l2x * l5x * vel[0] * sin06) - (
                              l2x * l5x * mg * vel[0] * sin06) - (c5x * c5z * l5m * vel[0] * cos06 / 2.0) - (
                              cgx * cgz * mg * vel[0] * cos06 / 2.0) + (c5x * l5m * l4z * vel[0] * cos07 / 4.0) + (
                              cgx * l4z * mg * vel[0] * cos07 / 4.0) - (c5y * l5m * l5x * vel[0] * cos06 / 2.0) - (
                              c5x * l5m * l3z * vel[0] * cos22 / 4.0) - (c5x * l5m * l4z * vel[0] * cos21 / 4.0) - (
                              c5x * l5m * l5z * vel[0] * cos06 / 2.0) - (cgy * l5x * mg * vel[0] * cos06 / 2.0) - (
                              cgx * l3z * mg * vel[0] * cos22 / 4.0) - (cgx * l4z * mg * vel[0] * cos21 / 4.0) - (
                              cgx * l5z * mg * vel[0] * cos06 / 2.0) + (c5y * c5z * l5m * vel[0] * sin13 / 2.0) + (
                              cgy * cgz * mg * vel[0] * sin13 / 2.0) - (c5y * l5m * l4z * vel[0] * sin14 / 4.0) - (
                              cgy * l4z * mg * vel[0] * sin14 / 4.0) - (c5x * l5m * l5x * vel[0] * sin13 / 2.0) - (
                              c5y * l5m * l3z * vel[0] * sin27 / 4.0) - (c5y * l5m * l4z * vel[0] * sin25 / 4.0) + (
                              c5y * l5m * l5z * vel[0] * sin13 / 2.0) - (cgx * l5x * mg * vel[0] * sin13 / 2.0) - (
                              cgy * l3z * mg * vel[0] * sin27 / 4.0) - (cgy * l4z * mg * vel[0] * sin25 / 4.0) + (
                              cgy * l5z * mg * vel[0] * sin13 / 2.0) + (c4x * l4m * l4z * vel[0] * cos03 / 2.0) - (
                              c4x * l4m * l3z * vel[0] * cos20 / 2.0) - (c4x * l4m * l4z * vel[0] * cos19 / 2.0) + (
                              c5x * l5m * l3z * vel[0] * cos13 / 4.0) + (cgx * l3z * mg * vel[0] * cos13 / 4.0) + (
                              l5m * l5x * l4z * vel[0] * cos03 / 2.0) + (l5x * l4z * mg * vel[0] * cos03 / 2.0) - (
                              l5m * l5x * l3z * vel[0] * cos20 / 2.0) - (l5m * l5x * l4z * vel[0] * cos19 / 2.0) - (
                              c4z * l4m * l4z * vel[0] * sin08 / 2.0) - (c5z * l5m * l4z * vel[0] * sin08 / 2.0) - (
                              l5x * l3z * mg * vel[0] * cos20 / 2.0) - (l5x * l4z * mg * vel[0] * cos19 / 2.0) - (
                              cgz * l4z * mg * vel[0] * sin08 / 2.0) + (c5x * c5z * l5m * vel[1] * cos01 / 2.0) - (
                              c5x * c5z * l5m * vel[2] * cos01 / 2.0) + (c5x * c5z * l5m * vel[3] * cos01 / 2.0) + (
                              cgx * cgz * mg * vel[1] * cos01 / 2.0) - (cgx * cgz * mg * vel[2] * cos01 / 2.0) + (
                              cgx * cgz * mg * vel[3] * cos01 / 2.0) + (c5y * l5m * l3z * vel[0] * sin17 / 4.0) + (
                              c4z * l4m * l3z * vel[0] * sin19 / 2.0) + (c4z * l4m * l4z * vel[0] * sin20 / 2.0) + (
                              c5z * l5m * l3z * vel[0] * sin19 / 2.0) + (c5z * l5m * l4z * vel[0] * sin20 / 2.0) + (
                              cgy * l3z * mg * vel[0] * sin17 / 4.0) + (cgz * l3z * mg * vel[0] * sin19 / 2.0) + (
                              cgz * l4z * mg * vel[0] * sin20 / 2.0) - (l5m * l4z * l5z * vel[0] * sin08 / 2.0)
        C[0, 4] = (I5x * vel[0] * sin09 / 8.0) - (I5x * vel[0] * sin10 / 8.0) + (I5y * vel[0] * sin10 / 8.0) - (
                    I5y * vel[0] * sin09 / 8.0) - (Igx * vel[0] * sin10 / 8.0) + (Igx * vel[0] * sin09 / 8.0) + (
                              Igy * vel[0] * sin10 / 8.0) - (Igy * vel[0] * sin09 / 8.0) - (
                              I5x * vel[0] * sin00 / 4.0) + (I5y * vel[0] * sin00 / 4.0) - (
                              Igx * vel[0] * sin00 / 4.0) + (Igy * vel[0] * sin00 / 4.0) - (
                              I5x * vel[1] * sin01 / 4.0) - (I5x * vel[1] * sin02 / 4.0) + (
                              I5x * vel[2] * sin01 / 4.0) + (I5x * vel[2] * sin02 / 4.0) - (
                              I5x * vel[3] * sin01 / 4.0) - (I5x * vel[3] * sin02 / 4.0) + (
                              I5y * vel[1] * sin01 / 4.0) + (I5y * vel[1] * sin02 / 4.0) - (
                              I5y * vel[2] * sin01 / 4.0) - (I5y * vel[2] * sin02 / 4.0) + (
                              I5y * vel[3] * sin01 / 4.0) + (I5y * vel[3] * sin02 / 4.0) - (
                              Igx * vel[1] * sin01 / 4.0) - (Igx * vel[1] * sin02 / 4.0) + (
                              Igx * vel[2] * sin01 / 4.0) + (Igx * vel[2] * sin02 / 4.0) - (
                              Igx * vel[3] * sin01 / 4.0) - (Igx * vel[3] * sin02 / 4.0) + (
                              Igy * vel[1] * sin01 / 4.0) + (Igy * vel[1] * sin02 / 4.0) - (
                              Igy * vel[2] * sin01 / 4.0) - (Igy * vel[2] * sin02 / 4.0) + (
                              Igy * vel[3] * sin01 / 4.0) + (Igy * vel[3] * sin02 / 4.0) - (
                              I5z * vel[1] * sin06 / 2.0) + (I5z * vel[2] * sin06 / 2.0) - (
                              I5z * vel[3] * sin06 / 2.0) - (Igz * vel[1] * sin06 / 2.0) + (
                              Igz * vel[2] * sin06 / 2.0) - (Igz * vel[3] * sin06 / 2.0) + (
                              pow(c5x, 2.0) * l5m * vel[0] * sin10 / 8.0) - (
                              pow(c5x, 2.0) * l5m * vel[0] * sin09 / 8.0) - (
                              pow(c5y, 2.0) * l5m * vel[0] * sin10 / 8.0) + (
                              pow(c5y, 2.0) * l5m * vel[0] * sin09 / 8.0) + (
                              pow(cgx, 2.0) * mg * vel[0] * sin10 / 8.0) - (
                              pow(cgx, 2.0) * mg * vel[0] * sin09 / 8.0) - (
                              pow(cgy, 2.0) * mg * vel[0] * sin10 / 8.0) + (
                              pow(cgy, 2.0) * mg * vel[0] * sin09 / 8.0) + (
                              pow(c5x, 2.0) * l5m * vel[0] * sin00 / 4.0) - (
                              pow(c5y, 2.0) * l5m * vel[0] * sin00 / 4.0) + (
                              pow(cgx, 2.0) * mg * vel[0] * sin00 / 4.0) - (
                              pow(cgy, 2.0) * mg * vel[0] * sin00 / 4.0) + (
                              pow(c5x, 2.0) * l5m * vel[1] * sin01 / 4.0) + (
                              pow(c5x, 2.0) * l5m * vel[1] * sin02 / 4.0) - (
                              pow(c5x, 2.0) * l5m * vel[2] * sin01 / 4.0) - (
                              pow(c5x, 2.0) * l5m * vel[2] * sin02 / 4.0) + (
                              pow(c5x, 2.0) * l5m * vel[3] * sin01 / 4.0) + (
                              pow(c5x, 2.0) * l5m * vel[3] * sin02 / 4.0) - (
                              pow(c5y, 2.0) * l5m * vel[1] * sin01 / 4.0) - (
                              pow(c5y, 2.0) * l5m * vel[1] * sin02 / 4.0) + (
                              pow(c5y, 2.0) * l5m * vel[2] * sin01 / 4.0) + (
                              pow(c5y, 2.0) * l5m * vel[2] * sin02 / 4.0) - (
                              pow(c5y, 2.0) * l5m * vel[3] * sin01 / 4.0) - (
                              pow(c5y, 2.0) * l5m * vel[3] * sin02 / 4.0) + (
                              pow(cgx, 2.0) * mg * vel[1] * sin01 / 4.0) + (
                              pow(cgx, 2.0) * mg * vel[1] * sin02 / 4.0) - (
                              pow(cgx, 2.0) * mg * vel[2] * sin01 / 4.0) - (
                              pow(cgx, 2.0) * mg * vel[2] * sin02 / 4.0) + (
                              pow(cgx, 2.0) * mg * vel[3] * sin01 / 4.0) + (
                              pow(cgx, 2.0) * mg * vel[3] * sin02 / 4.0) - (
                              pow(cgy, 2.0) * mg * vel[1] * sin01 / 4.0) - (
                              pow(cgy, 2.0) * mg * vel[1] * sin02 / 4.0) + (
                              pow(cgy, 2.0) * mg * vel[2] * sin01 / 4.0) + (
                              pow(cgy, 2.0) * mg * vel[2] * sin02 / 4.0) - (
                              pow(cgy, 2.0) * mg * vel[3] * sin01 / 4.0) - (
                              pow(cgy, 2.0) * mg * vel[3] * sin02 / 4.0) - (
                              pow(c5x, 2.0) * l5m * vel[1] * sin06 / 2.0) + (
                              pow(c5x, 2.0) * l5m * vel[2] * sin06 / 2.0) - (
                              pow(c5x, 2.0) * l5m * vel[3] * sin06 / 2.0) - (
                              pow(c5y, 2.0) * l5m * vel[1] * sin06 / 2.0) + (
                              pow(c5y, 2.0) * l5m * vel[2] * sin06 / 2.0) - (
                              pow(c5y, 2.0) * l5m * vel[3] * sin06 / 2.0) - (
                              pow(cgx, 2.0) * mg * vel[1] * sin06 / 2.0) + (
                              pow(cgx, 2.0) * mg * vel[2] * sin06 / 2.0) - (
                              pow(cgx, 2.0) * mg * vel[3] * sin06 / 2.0) - (
                              pow(cgy, 2.0) * mg * vel[1] * sin06 / 2.0) + (
                              pow(cgy, 2.0) * mg * vel[2] * sin06 / 2.0) - (
                              pow(cgy, 2.0) * mg * vel[3] * sin06 / 2.0) - (c5y * l5m * l2x * vel[0] * cos01 / 2.0) - (
                              c5y * l5m * l5x * vel[4] * cos01 / 2.0) - (c5x * l5m * l5z * vel[4] * cos01 / 2.0) - (
                              cgy * l2x * mg * vel[0] * cos01 / 2.0) - (cgy * l5x * mg * vel[4] * cos01 / 2.0) - (
                              cgx * l5z * mg * vel[4] * cos01 / 2.0) + (c5y * c5z * l5m * vel[4] * sin03 / 2.0) + (
                              cgy * cgz * mg * vel[4] * sin03 / 2.0) - (c5x * l5m * l2x * vel[0] * sin03 / 2.0) - (
                              c5x * l5m * l5x * vel[4] * sin03 / 2.0) + (c5y * l5m * l5z * vel[4] * sin03 / 2.0) - (
                              cgx * l2x * mg * vel[0] * sin03 / 2.0) - (cgx * l5x * mg * vel[4] * sin03 / 2.0) + (
                              cgy * l5z * mg * vel[4] * sin03 / 2.0) - (c5x * c5y * l5m * vel[0] * cos23 / 4.0) - (
                              c5x * c5y * l5m * vel[0] * cos26 / 4.0) + (c5x * c5z * l5m * vel[0] * cos05 / 4.0) - (
                              cgx * cgy * mg * vel[0] * cos23 / 4.0) - (cgx * cgy * mg * vel[0] * cos26 / 4.0) + (
                              cgx * cgz * mg * vel[0] * cos05 / 4.0) - (c5x * l5m * l4z * vel[0] * cos09 / 4.0) + (
                              c5x * l5m * l3z * vel[4] * cos29 / 2.0) - (cgx * l4z * mg * vel[0] * cos09 / 4.0) + (
                              cgx * l3z * mg * vel[4] * cos29 / 2.0) - (c5y * l5m * l5x * vel[0] * cos05 / 4.0) + (
                              c5x * l5m * l5z * vel[0] * cos05 / 4.0) - (cgy * l5x * mg * vel[0] * cos05 / 4.0) + (
                              cgx * l5z * mg * vel[0] * cos05 / 4.0) + (c5y * c5z * l5m * vel[0] * sin12 / 4.0) + (
                              cgy * cgz * mg * vel[0] * sin12 / 4.0) - (c5y * l5m * l4z * vel[0] * sin15 / 4.0) + (
                              c5y * l5m * l3z * vel[4] * sin32 / 2.0) - (cgy * l4z * mg * vel[0] * sin15 / 4.0) + (
                              cgy * l3z * mg * vel[4] * sin32 / 2.0) + (c5x * l5m * l5x * vel[0] * sin12 / 4.0) + (
                              c5y * l5m * l5z * vel[0] * sin12 / 4.0) + (cgx * l5x * mg * vel[0] * sin12 / 4.0) + (
                              cgy * l5z * mg * vel[0] * sin12 / 4.0) + (c5x * c5y * l5m * vel[0] * cos11 / 2.0) + (
                              cgx * cgy * mg * vel[0] * cos11 / 2.0) - (c5x * c5y * l5m * vel[1] * cos08 / 2.0) + (
                              c5x * c5y * l5m * vel[1] * cos10 / 2.0) + (c5x * c5y * l5m * vel[2] * cos08 / 2.0) - (
                              c5x * c5y * l5m * vel[2] * cos10 / 2.0) - (c5x * c5y * l5m * vel[3] * cos08 / 2.0) + (
                              c5x * c5y * l5m * vel[3] * cos10 / 2.0) + (c5x * c5z * l5m * vel[4] * cos00 / 2.0) - (
                              cgx * cgy * mg * vel[1] * cos08 / 2.0) + (cgx * cgy * mg * vel[1] * cos10 / 2.0) + (
                              cgx * cgy * mg * vel[2] * cos08 / 2.0) - (cgx * cgy * mg * vel[2] * cos10 / 2.0) - (
                              cgx * cgy * mg * vel[3] * cos08 / 2.0) + (cgx * cgy * mg * vel[3] * cos10 / 2.0) + (
                              cgx * cgz * mg * vel[4] * cos00 / 2.0) - (c5y * l5m * l2x * vel[0] * cos00 / 2.0) - (
                              c5y * l5m * l5x * vel[4] * cos00 / 2.0) - (c5x * l5m * l3z * vel[0] * cos25 / 4.0) - (
                              c5x * l5m * l4z * vel[0] * cos24 / 4.0) + (c5x * l5m * l5z * vel[4] * cos00 / 2.0) - (
                              cgy * l2x * mg * vel[0] * cos00 / 2.0) - (cgy * l5x * mg * vel[4] * cos00 / 2.0) - (
                              cgx * l3z * mg * vel[0] * cos25 / 4.0) - (cgx * l4z * mg * vel[0] * cos24 / 4.0) + (
                              cgx * l5z * mg * vel[4] * cos00 / 2.0) + (c5y * c5z * l5m * vel[4] * sin04 / 2.0) + (
                              cgy * cgz * mg * vel[4] * sin04 / 2.0) + (c5x * l5m * l2x * vel[0] * sin04 / 2.0) + (
                              c5x * l5m * l5x * vel[4] * sin04 / 2.0) + (c5y * l5m * l3z * vel[0] * sin26 / 4.0) + (
                              c5y * l5m * l4z * vel[0] * sin28 / 4.0) + (c5y * l5m * l5z * vel[4] * sin04 / 2.0) + (
                              cgx * l2x * mg * vel[0] * sin04 / 2.0) + (cgx * l5x * mg * vel[4] * sin04 / 2.0) + (
                              cgy * l3z * mg * vel[0] * sin26 / 4.0) + (cgy * l4z * mg * vel[0] * sin28 / 4.0) + (
                              cgy * l5z * mg * vel[4] * sin04 / 2.0) - (c5x * l5m * l3z * vel[0] * cos14 / 4.0) - (
                              c5x * l5m * l4z * vel[4] * cos16 / 2.0) - (cgx * l3z * mg * vel[0] * cos14 / 4.0) - (
                              cgx * l4z * mg * vel[4] * cos16 / 2.0) + (c5y * l5m * l3z * vel[0] * sin16 / 4.0) + (
                              c5y * l5m * l4z * vel[4] * sin21 / 2.0) + (cgy * l3z * mg * vel[0] * sin16 / 4.0) + (
                              cgy * l4z * mg * vel[4] * sin21 / 2.0) - (c5x * c5z * l5m * vel[0] * cos06 / 4.0) - (
                              cgx * cgz * mg * vel[0] * cos06 / 4.0) + (c5x * l5m * l4z * vel[0] * cos07 / 4.0) - (
                              c5x * l5m * l3z * vel[4] * cos28 / 2.0) + (cgx * l4z * mg * vel[0] * cos07 / 4.0) - (
                              cgx * l3z * mg * vel[4] * cos28 / 2.0) - (c5y * l5m * l5x * vel[0] * cos06 / 4.0) + (
                              c5x * l5m * l3z * vel[0] * cos22 / 4.0) + (c5x * l5m * l4z * vel[0] * cos21 / 4.0) - (
                              c5x * l5m * l5z * vel[0] * cos06 / 4.0) - (cgy * l5x * mg * vel[0] * cos06 / 4.0) + (
                              cgx * l3z * mg * vel[0] * cos22 / 4.0) + (cgx * l4z * mg * vel[0] * cos21 / 4.0) - (
                              cgx * l5z * mg * vel[0] * cos06 / 4.0) + (c5y * c5z * l5m * vel[0] * sin13 / 4.0) + (
                              cgy * cgz * mg * vel[0] * sin13 / 4.0) - (c5y * l5m * l4z * vel[0] * sin14 / 4.0) + (
                              c5y * l5m * l3z * vel[4] * sin31 / 2.0) - (cgy * l4z * mg * vel[0] * sin14 / 4.0) + (
                              cgy * l3z * mg * vel[4] * sin31 / 2.0) - (c5x * l5m * l5x * vel[0] * sin13 / 4.0) + (
                              c5y * l5m * l3z * vel[0] * sin27 / 4.0) + (c5y * l5m * l4z * vel[0] * sin25 / 4.0) + (
                              c5y * l5m * l5z * vel[0] * sin13 / 4.0) - (cgx * l5x * mg * vel[0] * sin13 / 4.0) + (
                              cgy * l3z * mg * vel[0] * sin27 / 4.0) + (cgy * l4z * mg * vel[0] * sin25 / 4.0) + (
                              cgy * l5z * mg * vel[0] * sin13 / 4.0) - (c5y * l5m * l5x * vel[0] * cos02 / 2.0) - (
                              c5y * l5m * l2x * vel[4] * cos02) - (cgy * l5x * mg * vel[0] * cos02 / 2.0) - (
                              cgy * l2x * mg * vel[4] * cos02) + (c5x * l5m * l3z * vel[0] * cos13 / 4.0) + (
                              c5x * l5m * l4z * vel[4] * cos17 / 2.0) + (cgx * l3z * mg * vel[0] * cos13 / 4.0) + (
                              cgx * l4z * mg * vel[4] * cos17 / 2.0) - (c5x * l5m * l5x * vel[0] * sin07 / 2.0) - (
                              c5x * l5m * l2x * vel[4] * sin07) - (cgx * l5x * mg * vel[0] * sin07 / 2.0) - (
                              cgx * l2x * mg * vel[4] * sin07) - (c5x * c5z * l5m * vel[4] * cos01 / 2.0) - (
                              cgx * cgz * mg * vel[4] * cos01 / 2.0) + (c5y * l5m * l3z * vel[0] * sin17 / 4.0) + (
                              c5y * l5m * l4z * vel[4] * sin22 / 2.0) + (cgy * l3z * mg * vel[0] * sin17 / 4.0) + (
                              cgy * l4z * mg * vel[4] * sin22 / 2.0)
        C[1, 0] = (I5y * vel[0] * sin10 / 8.0) - (I5x * vel[0] * sin09 / 8.0) - (I5x * vel[0] * sin10 / 8.0) + (
                    I5y * vel[0] * sin09 / 8.0) - (Igx * vel[0] * sin10 / 8.0) - (Igx * vel[0] * sin09 / 8.0) + (
                              Igy * vel[0] * sin10 / 8.0) + (Igy * vel[0] * sin09 / 8.0) - (
                              I2x * vel[0] * sin23 / 2.0) + (I2z * vel[0] * sin23 / 2.0) - (
                              I4x * vel[0] * sin05 / 2.0) - (I5x * vel[0] * sin05 / 4.0) - (
                              I5y * vel[0] * sin05 / 4.0) + (I4z * vel[0] * sin05 / 2.0) + (
                              I5z * vel[0] * sin05 / 2.0) - (Igx * vel[0] * sin05 / 4.0) - (
                              Igy * vel[0] * sin05 / 4.0) + (Igz * vel[0] * sin05 / 2.0) - (
                              I5x * vel[4] * sin01 / 4.0) - (I5x * vel[4] * sin02 / 4.0) + (
                              I5y * vel[4] * sin01 / 4.0) + (I5y * vel[4] * sin02 / 4.0) - (
                              Igx * vel[4] * sin01 / 4.0) - (Igx * vel[4] * sin02 / 4.0) + (
                              Igy * vel[4] * sin01 / 4.0) + (Igy * vel[4] * sin02 / 4.0) - (
                              I3x * vel[0] * sin18 / 2.0) + (I3z * vel[0] * sin18 / 2.0) + (
                              I5z * vel[4] * sin06 / 2.0) + (Igz * vel[4] * sin06 / 2.0) + (
                              pow(c5x, 2.0) * l5m * vel[0] * sin10 / 8.0) + (
                              pow(c5x, 2.0) * l5m * vel[0] * sin09 / 8.0) - (
                              pow(c5y, 2.0) * l5m * vel[0] * sin10 / 8.0) - (
                              pow(c5y, 2.0) * l5m * vel[0] * sin09 / 8.0) + (
                              pow(cgx, 2.0) * mg * vel[0] * sin10 / 8.0) + (
                              pow(cgx, 2.0) * mg * vel[0] * sin09 / 8.0) - (
                              pow(cgy, 2.0) * mg * vel[0] * sin10 / 8.0) - (
                              pow(cgy, 2.0) * mg * vel[0] * sin09 / 8.0) + (
                              pow(c2x, 2.0) * l2m * vel[0] * sin23 / 2.0) - (
                              pow(c2z, 2.0) * l2m * vel[0] * sin23 / 2.0) + (
                              pow(c4x, 2.0) * l4m * vel[0] * sin05 / 2.0) + (
                              pow(c5x, 2.0) * l5m * vel[0] * sin05 / 4.0) + (
                              pow(c5y, 2.0) * l5m * vel[0] * sin05 / 4.0) - (
                              pow(c4z, 2.0) * l4m * vel[0] * sin05 / 2.0) - (
                              pow(c5z, 2.0) * l5m * vel[0] * sin05 / 2.0) + (
                              pow(cgx, 2.0) * mg * vel[0] * sin05 / 4.0) + (
                              pow(cgy, 2.0) * mg * vel[0] * sin05 / 4.0) - (
                              pow(cgz, 2.0) * mg * vel[0] * sin05 / 2.0) - (
                              l3m * pow(l3z, 2.0) * vel[0] * sin23 / 2.0) - (
                              l4m * pow(l3z, 2.0) * vel[0] * sin23 / 2.0) - (
                              l5m * pow(l3z, 2.0) * vel[0] * sin23 / 2.0) - (
                              pow(l3z, 2.0) * mg * vel[0] * sin23 / 2.0) + (
                              l5m * pow(l5x, 2.0) * vel[0] * sin05 / 2.0) - (
                              l5m * pow(l5z, 2.0) * vel[0] * sin05 / 2.0) + (
                              pow(l5x, 2.0) * mg * vel[0] * sin05 / 2.0) - (
                              pow(l5z, 2.0) * mg * vel[0] * sin05 / 2.0) + (
                              pow(c5x, 2.0) * l5m * vel[4] * sin01 / 4.0) + (
                              pow(c5x, 2.0) * l5m * vel[4] * sin02 / 4.0) - (
                              pow(c5y, 2.0) * l5m * vel[4] * sin01 / 4.0) - (
                              pow(c5y, 2.0) * l5m * vel[4] * sin02 / 4.0) + (
                              pow(cgx, 2.0) * mg * vel[4] * sin01 / 4.0) + (
                              pow(cgx, 2.0) * mg * vel[4] * sin02 / 4.0) - (
                              pow(cgy, 2.0) * mg * vel[4] * sin01 / 4.0) - (
                              pow(cgy, 2.0) * mg * vel[4] * sin02 / 4.0) + (
                              pow(c3x, 2.0) * l3m * vel[0] * sin18 / 2.0) - (
                              pow(c3z, 2.0) * l3m * vel[0] * sin18 / 2.0) - (
                              l4m * pow(l4z, 2.0) * vel[0] * sin18 / 2.0) - (
                              l5m * pow(l4z, 2.0) * vel[0] * sin18 / 2.0) - (
                              pow(l4z, 2.0) * mg * vel[0] * sin18 / 2.0) + (
                              pow(c5x, 2.0) * l5m * vel[4] * sin06 / 2.0) + (
                              pow(c5y, 2.0) * l5m * vel[4] * sin06 / 2.0) + (
                              pow(cgx, 2.0) * mg * vel[4] * sin06 / 2.0) + (
                              pow(cgy, 2.0) * mg * vel[4] * sin06 / 2.0) + (c5y * l5m * l2x * vel[0] * cos01 / 2.0) + (
                              c5y * l5m * l5x * vel[4] * cos01 / 2.0) + (c5x * l5m * l5z * vel[4] * cos01 / 2.0) + (
                              cgy * l2x * mg * vel[0] * cos01 / 2.0) + (cgy * l5x * mg * vel[4] * cos01 / 2.0) + (
                              cgx * l5z * mg * vel[4] * cos01 / 2.0) - (l5m * l3z * l5z * vel[0] * sin19) - (
                              l5m * l4z * l5z * vel[0] * sin20) - (l3z * l5z * mg * vel[0] * sin19) - (
                              l4z * l5z * mg * vel[0] * sin20) - (c5y * c5z * l5m * vel[4] * sin03 / 2.0) - (
                              cgy * cgz * mg * vel[4] * sin03 / 2.0) + (c5x * l5m * l2x * vel[0] * sin03 / 2.0) + (
                              c5x * l5m * l5x * vel[4] * sin03 / 2.0) - (c5y * l5m * l5z * vel[4] * sin03 / 2.0) + (
                              cgx * l2x * mg * vel[0] * sin03 / 2.0) + (cgx * l5x * mg * vel[4] * sin03 / 2.0) - (
                              cgy * l5z * mg * vel[4] * sin03 / 2.0) - (c5x * c5y * l5m * vel[0] * cos23 / 4.0) + (
                              c5x * c5y * l5m * vel[0] * cos26 / 4.0) + (c5x * c5z * l5m * vel[0] * cos05 / 2.0) - (
                              cgx * cgy * mg * vel[0] * cos23 / 4.0) + (cgx * cgy * mg * vel[0] * cos26 / 4.0) + (
                              cgx * cgz * mg * vel[0] * cos05 / 2.0) + (c3z * l3m * l2x * vel[0] * cos18) + (
                              c5x * l5m * l3z * vel[4] * cos29 / 2.0) + (cgx * l3z * mg * vel[4] * cos29 / 2.0) - (
                              c5y * l5m * l5x * vel[0] * cos05 / 2.0) + (c5x * l5m * l5z * vel[0] * cos05 / 2.0) - (
                              cgy * l5x * mg * vel[0] * cos05 / 2.0) + (cgx * l5z * mg * vel[0] * cos05 / 2.0) + (
                              l4m * l2x * l4z * vel[0] * cos18) + (l5m * l2x * l4z * vel[0] * cos18) + (
                              l2x * l4z * mg * vel[0] * cos18) + (c5y * c5z * l5m * vel[0] * sin12 / 2.0) + (
                              cgy * cgz * mg * vel[0] * sin12 / 2.0) + (c3x * l3m * l2x * vel[0] * sin30) + (
                              c5y * l5m * l3z * vel[4] * sin32 / 2.0) + (cgy * l3z * mg * vel[4] * sin32 / 2.0) + (
                              c5x * l5m * l5x * vel[0] * sin12 / 2.0) + (c5y * l5m * l5z * vel[0] * sin12 / 2.0) + (
                              cgx * l5x * mg * vel[0] * sin12 / 2.0) + (cgy * l5z * mg * vel[0] * sin12 / 2.0) + (
                              c2x * c2z * l2m * vel[0] * cos33) + (c4x * c4z * l4m * vel[0] * cos15) + (
                              c5z * l5m * l5x * vel[0] * cos15) + (cgz * l5x * mg * vel[0] * cos15) + (
                              l5m * l5x * l5z * vel[0] * cos15) + (l5x * l5z * mg * vel[0] * cos15) - (
                              c5x * c5y * l5m * vel[4] * cos08 / 2.0) + (c5x * c5y * l5m * vel[4] * cos10 / 2.0) + (
                              c5x * c5z * l5m * vel[4] * cos00 / 2.0) - (cgx * cgy * mg * vel[4] * cos08 / 2.0) + (
                              cgx * cgy * mg * vel[4] * cos10 / 2.0) + (cgx * cgz * mg * vel[4] * cos00 / 2.0) - (
                              c5z * l5m * l5z * vel[0] * sin05) - (cgz * l5z * mg * vel[0] * sin05) - (
                              c5y * l5m * l2x * vel[0] * cos00 / 2.0) - (c5y * l5m * l5x * vel[4] * cos00 / 2.0) + (
                              c5x * l5m * l3z * vel[0] * cos25 / 2.0) + (c5x * l5m * l4z * vel[0] * cos24 / 2.0) + (
                              c5x * l5m * l5z * vel[4] * cos00 / 2.0) - (cgy * l2x * mg * vel[0] * cos00 / 2.0) - (
                              cgy * l5x * mg * vel[4] * cos00 / 2.0) + (cgx * l3z * mg * vel[0] * cos25 / 2.0) + (
                              cgx * l4z * mg * vel[0] * cos24 / 2.0) + (cgx * l5z * mg * vel[4] * cos00 / 2.0) + (
                              c5y * c5z * l5m * vel[4] * sin04 / 2.0) + (cgy * cgz * mg * vel[4] * sin04 / 2.0) + (
                              c3x * c3z * l3m * vel[0] * cos31) + (c5x * l5m * l2x * vel[0] * sin04 / 2.0) + (
                              c5x * l5m * l5x * vel[4] * sin04 / 2.0) - (c5y * l5m * l3z * vel[0] * sin26 / 2.0) - (
                              c5y * l5m * l4z * vel[0] * sin28 / 2.0) + (c5y * l5m * l5z * vel[4] * sin04 / 2.0) + (
                              cgx * l2x * mg * vel[0] * sin04 / 2.0) + (cgx * l5x * mg * vel[4] * sin04 / 2.0) - (
                              cgy * l3z * mg * vel[0] * sin26 / 2.0) - (cgy * l4z * mg * vel[0] * sin28 / 2.0) + (
                              cgy * l5z * mg * vel[4] * sin04 / 2.0) + (c3x * l3m * l3z * vel[0] * cos32) - (
                              c3z * l3m * l3z * vel[0] * sin29) + (c4z * l4m * l2x * vel[0] * cos04) + (
                              c5z * l5m * l2x * vel[0] * cos04) + (c5x * l5m * l4z * vel[4] * cos16 / 2.0) + (
                              cgz * l2x * mg * vel[0] * cos04) + (cgx * l4z * mg * vel[4] * cos16 / 2.0) - (
                              l4m * l3z * l4z * vel[0] * sin29) - (l5m * l3z * l4z * vel[0] * sin29) - (
                              l3z * l4z * mg * vel[0] * sin29) + (l5m * l2x * l5z * vel[0] * cos04) + (
                              l2x * l5z * mg * vel[0] * cos04) + (c4x * l4m * l2x * vel[0] * sin06) - (
                              c5y * l5m * l4z * vel[4] * sin21 / 2.0) - (cgy * l4z * mg * vel[4] * sin21 / 2.0) + (
                              l5m * l2x * l5x * vel[0] * sin06) + (l2x * l5x * mg * vel[0] * sin06) + (
                              c5x * c5z * l5m * vel[0] * cos06 / 2.0) + (cgx * cgz * mg * vel[0] * cos06 / 2.0) + (
                              c5x * l5m * l3z * vel[4] * cos28 / 2.0) + (cgx * l3z * mg * vel[4] * cos28 / 2.0) + (
                              c5y * l5m * l5x * vel[0] * cos06 / 2.0) + (c5x * l5m * l3z * vel[0] * cos22 / 2.0) + (
                              c5x * l5m * l4z * vel[0] * cos21 / 2.0) + (c5x * l5m * l5z * vel[0] * cos06 / 2.0) + (
                              cgy * l5x * mg * vel[0] * cos06 / 2.0) + (cgx * l3z * mg * vel[0] * cos22 / 2.0) + (
                              cgx * l4z * mg * vel[0] * cos21 / 2.0) + (cgx * l5z * mg * vel[0] * cos06 / 2.0) - (
                              c5y * c5z * l5m * vel[0] * sin13 / 2.0) - (cgy * cgz * mg * vel[0] * sin13 / 2.0) - (
                              c5y * l5m * l3z * vel[4] * sin31 / 2.0) - (cgy * l3z * mg * vel[4] * sin31 / 2.0) + (
                              c5x * l5m * l5x * vel[0] * sin13 / 2.0) + (c5y * l5m * l3z * vel[0] * sin27 / 2.0) + (
                              c5y * l5m * l4z * vel[0] * sin25 / 2.0) - (c5y * l5m * l5z * vel[0] * sin13 / 2.0) + (
                              cgx * l5x * mg * vel[0] * sin13 / 2.0) + (cgy * l3z * mg * vel[0] * sin27 / 2.0) + (
                              cgy * l4z * mg * vel[0] * sin25 / 2.0) - (cgy * l5z * mg * vel[0] * sin13 / 2.0) + (
                              c2z * l2m * l2x * vel[0] * cos27) + (c4x * l4m * l3z * vel[0] * cos20) + (
                              c4x * l4m * l4z * vel[0] * cos19) + (c5x * l5m * l4z * vel[4] * cos17 / 2.0) + (
                              cgx * l4z * mg * vel[4] * cos17 / 2.0) + (l3m * l2x * l3z * vel[0] * cos27) + (
                              l4m * l2x * l3z * vel[0] * cos27) + (l5m * l2x * l3z * vel[0] * cos27) + (
                              l2x * l3z * mg * vel[0] * cos27) + (l5m * l5x * l3z * vel[0] * cos20) + (
                              l5m * l5x * l4z * vel[0] * cos19) + (c2x * l2m * l2x * vel[0] * sin33) + (
                              l5x * l3z * mg * vel[0] * cos20) + (l5x * l4z * mg * vel[0] * cos19) + (
                              c5x * c5z * l5m * vel[4] * cos01 / 2.0) + (cgx * cgz * mg * vel[4] * cos01 / 2.0) + (
                              c5y * l5m * l4z * vel[4] * sin22 / 2.0) - (c4z * l4m * l3z * vel[0] * sin19) - (
                              c4z * l4m * l4z * vel[0] * sin20) - (c5z * l5m * l3z * vel[0] * sin19) - (
                              c5z * l5m * l4z * vel[0] * sin20) + (cgy * l4z * mg * vel[4] * sin22 / 2.0) - (
                              cgz * l3z * mg * vel[0] * sin19) - (cgz * l4z * mg * vel[0] * sin20)
        C[1, 1] = -(vel[4] * ((I5y * sin00) - (I5x * sin00) - (Igx * sin00) + (Igy * sin00) + (
                    2.0 * l5m * ((c5y * cos02) + (c5x * sin07)) * (
                        l5x + (c5x * cos02) - (c5y * sin07) + (l4z * sin08) - (l3z * sin11))) + (
                                          2.0 * mg * ((cgy * cos02) + (cgx * sin07)) * (
                                              l5x + (cgx * cos02) - (cgy * sin07) + (l4z * sin08) - (
                                                  l3z * sin11)))) / 2.0) - (vel[3] * (
                    (2.0 * l4m * ((l4z * sin08) - (l3z * sin11)) * (c4z + (l4z * cos03) + (l3z * cos12))) - (
                        2.0 * mg * ((l4z * cos03) + (l3z * cos12)) * (
                            l5x + (cgx * cos02) - (cgy * sin07) + (l4z * sin08) - (l3z * sin11))) - (
                                2.0 * l4m * ((l4z * cos03) + (l3z * cos12)) * (c4x + (l4z * sin08) - (l3z * sin11))) - (
                                2.0 * l5m * ((l4z * cos03) + (l3z * cos12)) * (
                                    l5x + (c5x * cos02) - (c5y * sin07) + (l4z * sin08) - (l3z * sin11))) + (
                                2.0 * l5m * pow(cos02, 2.0) * ((l4z * sin08) - (l3z * sin11)) * (
                                    c5z + l5z + (l4z * cos03) + (l3z * cos12))) + (
                                2.0 * mg * pow(cos02, 2.0) * ((l4z * sin08) - (l3z * sin11)) * (
                                    cgz + l5z + (l4z * cos03) + (l3z * cos12))) + (
                                2.0 * l5m * pow(sin07, 2.0) * ((l4z * sin08) - (l3z * sin11)) * (
                                    c5z + l5z + (l4z * cos03) + (l3z * cos12))) + (
                                2.0 * mg * pow(sin07, 2.0) * ((l4z * sin08) - (l3z * sin11)) * (
                                    cgz + l5z + (l4z * cos03) + (l3z * cos12)))) / 2.0) - (l3z * vel[2] * (
                    (2.0 * c3x * l3m * cos30) + (c5x * l5m * cos13) + (cgx * mg * cos13) + (2.0 * c3z * l3m * sin24) + (
                        c5y * l5m * sin17) + (cgy * mg * sin17) + (2.0 * l4m * l4z * sin24) + (
                                2.0 * l5m * l4z * sin24) + (2.0 * l4z * mg * sin24) + (2.0 * c4x * l4m * cos12) + (
                                2.0 * l5m * l5x * cos12) + (2.0 * l5x * mg * cos12) + (2.0 * c4z * l4m * sin11) + (
                                2.0 * c5z * l5m * sin11) + (2.0 * cgz * mg * sin11) + (2.0 * l5m * l5z * sin11) + (
                                2.0 * l5z * mg * sin11) + (c5x * l5m * cos14) + (cgx * mg * cos14) - (
                                c5y * l5m * sin16) - (cgy * mg * sin16)) / 2.0)
        C[1, 2] = (I5y * vel[4] * sin00 / 2.0) - (I5x * vel[4] * sin00 / 2.0) - (Igx * vel[4] * sin00 / 2.0) + (
                    Igy * vel[4] * sin00 / 2.0) + (pow(c5x, 2.0) * l5m * vel[4] * sin00 / 2.0) - (
                              pow(c5y, 2.0) * l5m * vel[4] * sin00 / 2.0) + (
                              pow(cgx, 2.0) * mg * vel[4] * sin00 / 2.0) - (
                              pow(cgy, 2.0) * mg * vel[4] * sin00 / 2.0) - (l3z * l4z * mg * vel[1] * sin24) + (
                              l3z * l4z * mg * vel[2] * sin24) + (l4z * l5z * mg * vel[3] * sin08) - (
                              c4x * l4m * l3z * vel[1] * cos12) + (c4x * l4m * l3z * vel[2] * cos12) - (
                              c4x * l4m * l3z * vel[3] * cos12) - (c5x * l5m * l4z * vel[3] * cos09 / 2.0) + (
                              c5x * l5m * l4z * vel[4] * cos09 / 2.0) - (cgx * l4z * mg * vel[3] * cos09 / 2.0) + (
                              cgx * l4z * mg * vel[4] * cos09 / 2.0) - (l5m * l5x * l3z * vel[1] * cos12) + (
                              l5m * l5x * l3z * vel[2] * cos12) - (l5m * l5x * l3z * vel[3] * cos12) - (
                              l5x * l3z * mg * vel[1] * cos12) + (l5x * l3z * mg * vel[2] * cos12) - (
                              l5x * l3z * mg * vel[3] * cos12) - (c5y * l5m * l4z * vel[3] * sin15 / 2.0) + (
                              c5y * l5m * l4z * vel[4] * sin15 / 2.0) - (c4z * l4m * l3z * vel[1] * sin11) + (
                              c4z * l4m * l3z * vel[2] * sin11) - (c4z * l4m * l3z * vel[3] * sin11) - (
                              c5z * l5m * l3z * vel[1] * sin11) + (c5z * l5m * l3z * vel[2] * sin11) - (
                              c5z * l5m * l3z * vel[3] * sin11) - (cgy * l4z * mg * vel[3] * sin15 / 2.0) + (
                              cgy * l4z * mg * vel[4] * sin15 / 2.0) - (cgz * l3z * mg * vel[1] * sin11) + (
                              cgz * l3z * mg * vel[2] * sin11) - (cgz * l3z * mg * vel[3] * sin11) - (
                              l5m * l3z * l5z * vel[1] * sin11) + (l5m * l3z * l5z * vel[2] * sin11) - (
                              l5m * l3z * l5z * vel[3] * sin11) - (l3z * l5z * mg * vel[1] * sin11) + (
                              l3z * l5z * mg * vel[2] * sin11) - (l3z * l5z * mg * vel[3] * sin11) + (
                              c5x * c5y * l5m * vel[4] * cos11) + (cgx * cgy * mg * vel[4] * cos11) - (
                              c5x * l5m * l3z * vel[1] * cos14 / 2.0) + (c5x * l5m * l3z * vel[2] * cos14 / 2.0) - (
                              c5x * l5m * l3z * vel[3] * cos14 / 2.0) + (c5x * l5m * l3z * vel[4] * cos14 / 2.0) - (
                              cgx * l3z * mg * vel[1] * cos14 / 2.0) + (cgx * l3z * mg * vel[2] * cos14 / 2.0) - (
                              cgx * l3z * mg * vel[3] * cos14 / 2.0) + (cgx * l3z * mg * vel[4] * cos14 / 2.0) + (
                              c5y * l5m * l3z * vel[1] * sin16 / 2.0) - (c5y * l5m * l3z * vel[2] * sin16 / 2.0) + (
                              c5y * l5m * l3z * vel[3] * sin16 / 2.0) - (c5y * l5m * l3z * vel[4] * sin16 / 2.0) + (
                              cgy * l3z * mg * vel[1] * sin16 / 2.0) - (cgy * l3z * mg * vel[2] * sin16 / 2.0) + (
                              cgy * l3z * mg * vel[3] * sin16 / 2.0) - (cgy * l3z * mg * vel[4] * sin16 / 2.0) - (
                              c5x * l5m * l4z * vel[3] * cos07 / 2.0) - (c5x * l5m * l4z * vel[4] * cos07 / 2.0) - (
                              cgx * l4z * mg * vel[3] * cos07 / 2.0) - (cgx * l4z * mg * vel[4] * cos07 / 2.0) + (
                              c5y * l5m * l4z * vel[3] * sin14 / 2.0) + (c5y * l5m * l4z * vel[4] * sin14 / 2.0) + (
                              cgy * l4z * mg * vel[3] * sin14 / 2.0) + (cgy * l4z * mg * vel[4] * sin14 / 2.0) + (
                              c5y * l5m * l5x * vel[4] * cos02) - (c3x * l3m * l3z * vel[1] * cos30) + (
                              c3x * l3m * l3z * vel[2] * cos30) - (c4x * l4m * l4z * vel[3] * cos03) + (
                              cgy * l5x * mg * vel[4] * cos02) - (c5x * l5m * l3z * vel[1] * cos13 / 2.0) + (
                              c5x * l5m * l3z * vel[2] * cos13 / 2.0) - (c5x * l5m * l3z * vel[3] * cos13 / 2.0) - (
                              c5x * l5m * l3z * vel[4] * cos13 / 2.0) - (cgx * l3z * mg * vel[1] * cos13 / 2.0) + (
                              cgx * l3z * mg * vel[2] * cos13 / 2.0) - (cgx * l3z * mg * vel[3] * cos13 / 2.0) - (
                              cgx * l3z * mg * vel[4] * cos13 / 2.0) - (l5m * l5x * l4z * vel[3] * cos03) - (
                              l5x * l4z * mg * vel[3] * cos03) + (c5x * l5m * l5x * vel[4] * sin07) - (
                              c3z * l3m * l3z * vel[1] * sin24) + (c3z * l3m * l3z * vel[2] * sin24) + (
                              c4z * l4m * l4z * vel[3] * sin08) + (c5z * l5m * l4z * vel[3] * sin08) + (
                              cgx * l5x * mg * vel[4] * sin07) + (cgz * l4z * mg * vel[3] * sin08) - (
                              c5y * l5m * l3z * vel[1] * sin17 / 2.0) + (c5y * l5m * l3z * vel[2] * sin17 / 2.0) - (
                              c5y * l5m * l3z * vel[3] * sin17 / 2.0) - (c5y * l5m * l3z * vel[4] * sin17 / 2.0) - (
                              cgy * l3z * mg * vel[1] * sin17 / 2.0) + (cgy * l3z * mg * vel[2] * sin17 / 2.0) - (
                              cgy * l3z * mg * vel[3] * sin17 / 2.0) - (cgy * l3z * mg * vel[4] * sin17 / 2.0) - (
                              l4m * l3z * l4z * vel[1] * sin24) + (l4m * l3z * l4z * vel[2] * sin24) - (
                              l5m * l3z * l4z * vel[1] * sin24) + (l5m * l3z * l4z * vel[2] * sin24) + (
                              l5m * l4z * l5z * vel[3] * sin08)
        C[1, 3] = (I5x * vel[4] * sin00 / 2.0) - (I5y * vel[4] * sin00 / 2.0) + (Igx * vel[4] * sin00 / 2.0) - (
                    Igy * vel[4] * sin00 / 2.0) - (pow(c5x, 2.0) * l5m * vel[4] * sin00 / 2.0) + (
                              pow(c5y, 2.0) * l5m * vel[4] * sin00 / 2.0) - (
                              pow(cgx, 2.0) * mg * vel[4] * sin00 / 2.0) + (
                              pow(cgy, 2.0) * mg * vel[4] * sin00 / 2.0) - (l4z * l5z * mg * vel[1] * sin08) + (
                              l4z * l5z * mg * vel[2] * sin08) - (l4z * l5z * mg * vel[3] * sin08) + (
                              c4x * l4m * l3z * vel[1] * cos12) - (c4x * l4m * l3z * vel[2] * cos12) + (
                              c4x * l4m * l3z * vel[3] * cos12) + (c5x * l5m * l4z * vel[1] * cos09 / 2.0) - (
                              c5x * l5m * l4z * vel[2] * cos09 / 2.0) + (c5x * l5m * l4z * vel[3] * cos09 / 2.0) - (
                              c5x * l5m * l4z * vel[4] * cos09 / 2.0) + (cgx * l4z * mg * vel[1] * cos09 / 2.0) - (
                              cgx * l4z * mg * vel[2] * cos09 / 2.0) + (cgx * l4z * mg * vel[3] * cos09 / 2.0) - (
                              cgx * l4z * mg * vel[4] * cos09 / 2.0) + (l5m * l5x * l3z * vel[1] * cos12) - (
                              l5m * l5x * l3z * vel[2] * cos12) + (l5m * l5x * l3z * vel[3] * cos12) + (
                              l5x * l3z * mg * vel[1] * cos12) - (l5x * l3z * mg * vel[2] * cos12) + (
                              l5x * l3z * mg * vel[3] * cos12) + (c5y * l5m * l4z * vel[1] * sin15 / 2.0) - (
                              c5y * l5m * l4z * vel[2] * sin15 / 2.0) + (c5y * l5m * l4z * vel[3] * sin15 / 2.0) - (
                              c5y * l5m * l4z * vel[4] * sin15 / 2.0) + (c4z * l4m * l3z * vel[1] * sin11) - (
                              c4z * l4m * l3z * vel[2] * sin11) + (c4z * l4m * l3z * vel[3] * sin11) + (
                              c5z * l5m * l3z * vel[1] * sin11) - (c5z * l5m * l3z * vel[2] * sin11) + (
                              c5z * l5m * l3z * vel[3] * sin11) + (cgy * l4z * mg * vel[1] * sin15 / 2.0) - (
                              cgy * l4z * mg * vel[2] * sin15 / 2.0) + (cgy * l4z * mg * vel[3] * sin15 / 2.0) - (
                              cgy * l4z * mg * vel[4] * sin15 / 2.0) + (cgz * l3z * mg * vel[1] * sin11) - (
                              cgz * l3z * mg * vel[2] * sin11) + (cgz * l3z * mg * vel[3] * sin11) + (
                              l5m * l3z * l5z * vel[1] * sin11) - (l5m * l3z * l5z * vel[2] * sin11) + (
                              l5m * l3z * l5z * vel[3] * sin11) + (l3z * l5z * mg * vel[1] * sin11) - (
                              l3z * l5z * mg * vel[2] * sin11) + (l3z * l5z * mg * vel[3] * sin11) - (
                              c5x * c5y * l5m * vel[4] * cos11) - (cgx * cgy * mg * vel[4] * cos11) + (
                              c5x * l5m * l3z * vel[1] * cos14 / 2.0) - (c5x * l5m * l3z * vel[2] * cos14 / 2.0) + (
                              c5x * l5m * l3z * vel[3] * cos14 / 2.0) - (c5x * l5m * l3z * vel[4] * cos14 / 2.0) + (
                              cgx * l3z * mg * vel[1] * cos14 / 2.0) - (cgx * l3z * mg * vel[2] * cos14 / 2.0) + (
                              cgx * l3z * mg * vel[3] * cos14 / 2.0) - (cgx * l3z * mg * vel[4] * cos14 / 2.0) - (
                              c5y * l5m * l3z * vel[1] * sin16 / 2.0) + (c5y * l5m * l3z * vel[2] * sin16 / 2.0) - (
                              c5y * l5m * l3z * vel[3] * sin16 / 2.0) + (c5y * l5m * l3z * vel[4] * sin16 / 2.0) - (
                              cgy * l3z * mg * vel[1] * sin16 / 2.0) + (cgy * l3z * mg * vel[2] * sin16 / 2.0) - (
                              cgy * l3z * mg * vel[3] * sin16 / 2.0) + (cgy * l3z * mg * vel[4] * sin16 / 2.0) + (
                              c5x * l5m * l4z * vel[1] * cos07 / 2.0) - (c5x * l5m * l4z * vel[2] * cos07 / 2.0) + (
                              c5x * l5m * l4z * vel[3] * cos07 / 2.0) + (c5x * l5m * l4z * vel[4] * cos07 / 2.0) + (
                              cgx * l4z * mg * vel[1] * cos07 / 2.0) - (cgx * l4z * mg * vel[2] * cos07 / 2.0) + (
                              cgx * l4z * mg * vel[3] * cos07 / 2.0) + (cgx * l4z * mg * vel[4] * cos07 / 2.0) - (
                              c5y * l5m * l4z * vel[1] * sin14 / 2.0) + (c5y * l5m * l4z * vel[2] * sin14 / 2.0) - (
                              c5y * l5m * l4z * vel[3] * sin14 / 2.0) - (c5y * l5m * l4z * vel[4] * sin14 / 2.0) - (
                              cgy * l4z * mg * vel[1] * sin14 / 2.0) + (cgy * l4z * mg * vel[2] * sin14 / 2.0) - (
                              cgy * l4z * mg * vel[3] * sin14 / 2.0) - (cgy * l4z * mg * vel[4] * sin14 / 2.0) - (
                              c5y * l5m * l5x * vel[4] * cos02) + (c4x * l4m * l4z * vel[1] * cos03) - (
                              c4x * l4m * l4z * vel[2] * cos03) + (c4x * l4m * l4z * vel[3] * cos03) - (
                              cgy * l5x * mg * vel[4] * cos02) + (c5x * l5m * l3z * vel[1] * cos13 / 2.0) - (
                              c5x * l5m * l3z * vel[2] * cos13 / 2.0) + (c5x * l5m * l3z * vel[3] * cos13 / 2.0) + (
                              c5x * l5m * l3z * vel[4] * cos13 / 2.0) + (cgx * l3z * mg * vel[1] * cos13 / 2.0) - (
                              cgx * l3z * mg * vel[2] * cos13 / 2.0) + (cgx * l3z * mg * vel[3] * cos13 / 2.0) + (
                              cgx * l3z * mg * vel[4] * cos13 / 2.0) + (l5m * l5x * l4z * vel[1] * cos03) - (
                              l5m * l5x * l4z * vel[2] * cos03) + (l5m * l5x * l4z * vel[3] * cos03) + (
                              l5x * l4z * mg * vel[1] * cos03) - (l5x * l4z * mg * vel[2] * cos03) + (
                              l5x * l4z * mg * vel[3] * cos03) - (c5x * l5m * l5x * vel[4] * sin07) - (
                              c4z * l4m * l4z * vel[1] * sin08) + (c4z * l4m * l4z * vel[2] * sin08) - (
                              c4z * l4m * l4z * vel[3] * sin08) - (c5z * l5m * l4z * vel[1] * sin08) + (
                              c5z * l5m * l4z * vel[2] * sin08) - (c5z * l5m * l4z * vel[3] * sin08) - (
                              cgx * l5x * mg * vel[4] * sin07) - (cgz * l4z * mg * vel[1] * sin08) + (
                              cgz * l4z * mg * vel[2] * sin08) - (cgz * l4z * mg * vel[3] * sin08) + (
                              c5y * l5m * l3z * vel[1] * sin17 / 2.0) - (c5y * l5m * l3z * vel[2] * sin17 / 2.0) + (
                              c5y * l5m * l3z * vel[3] * sin17 / 2.0) + (c5y * l5m * l3z * vel[4] * sin17 / 2.0) + (
                              cgy * l3z * mg * vel[1] * sin17 / 2.0) - (cgy * l3z * mg * vel[2] * sin17 / 2.0) + (
                              cgy * l3z * mg * vel[3] * sin17 / 2.0) + (cgy * l3z * mg * vel[4] * sin17 / 2.0) - (
                              l5m * l4z * l5z * vel[1] * sin08) + (l5m * l4z * l5z * vel[2] * sin08) - (
                              l5m * l4z * l5z * vel[3] * sin08)
        C[1, 4] = (I5x * vel[1] * sin00 / 2.0) - (I5x * vel[2] * sin00 / 2.0) + (I5x * vel[3] * sin00 / 2.0) - (
                    I5y * vel[1] * sin00 / 2.0) + (I5y * vel[2] * sin00 / 2.0) - (I5y * vel[3] * sin00 / 2.0) + (
                              Igx * vel[1] * sin00 / 2.0) - (Igx * vel[2] * sin00 / 2.0) + (
                              Igx * vel[3] * sin00 / 2.0) - (Igy * vel[1] * sin00 / 2.0) + (
                              Igy * vel[2] * sin00 / 2.0) - (Igy * vel[3] * sin00 / 2.0) - (
                              I5x * vel[0] * sin01 / 4.0) - (I5x * vel[0] * sin02 / 4.0) + (
                              I5y * vel[0] * sin01 / 4.0) + (I5y * vel[0] * sin02 / 4.0) - (
                              Igx * vel[0] * sin01 / 4.0) - (Igx * vel[0] * sin02 / 4.0) + (
                              Igy * vel[0] * sin01 / 4.0) + (Igy * vel[0] * sin02 / 4.0) + (
                              I5z * vel[0] * sin06 / 2.0) + (Igz * vel[0] * sin06 / 2.0) - (
                              pow(c5x, 2.0) * l5m * vel[1] * sin00 / 2.0) + (
                              pow(c5x, 2.0) * l5m * vel[2] * sin00 / 2.0) - (
                              pow(c5x, 2.0) * l5m * vel[3] * sin00 / 2.0) + (
                              pow(c5y, 2.0) * l5m * vel[1] * sin00 / 2.0) - (
                              pow(c5y, 2.0) * l5m * vel[2] * sin00 / 2.0) + (
                              pow(c5y, 2.0) * l5m * vel[3] * sin00 / 2.0) - (
                              pow(cgx, 2.0) * mg * vel[1] * sin00 / 2.0) + (
                              pow(cgx, 2.0) * mg * vel[2] * sin00 / 2.0) - (
                              pow(cgx, 2.0) * mg * vel[3] * sin00 / 2.0) + (
                              pow(cgy, 2.0) * mg * vel[1] * sin00 / 2.0) - (
                              pow(cgy, 2.0) * mg * vel[2] * sin00 / 2.0) + (
                              pow(cgy, 2.0) * mg * vel[3] * sin00 / 2.0) + (
                              pow(c5x, 2.0) * l5m * vel[0] * sin01 / 4.0) + (
                              pow(c5x, 2.0) * l5m * vel[0] * sin02 / 4.0) - (
                              pow(c5y, 2.0) * l5m * vel[0] * sin01 / 4.0) - (
                              pow(c5y, 2.0) * l5m * vel[0] * sin02 / 4.0) + (
                              pow(cgx, 2.0) * mg * vel[0] * sin01 / 4.0) + (
                              pow(cgx, 2.0) * mg * vel[0] * sin02 / 4.0) - (
                              pow(cgy, 2.0) * mg * vel[0] * sin01 / 4.0) - (
                              pow(cgy, 2.0) * mg * vel[0] * sin02 / 4.0) + (
                              pow(c5x, 2.0) * l5m * vel[0] * sin06 / 2.0) + (
                              pow(c5y, 2.0) * l5m * vel[0] * sin06 / 2.0) + (
                              pow(cgx, 2.0) * mg * vel[0] * sin06 / 2.0) + (
                              pow(cgy, 2.0) * mg * vel[0] * sin06 / 2.0) + (c5y * l5m * l5x * vel[0] * cos01 / 2.0) + (
                              c5x * l5m * l5z * vel[0] * cos01 / 2.0) + (cgy * l5x * mg * vel[0] * cos01 / 2.0) + (
                              cgx * l5z * mg * vel[0] * cos01 / 2.0) - (c5y * c5z * l5m * vel[0] * sin03 / 2.0) - (
                              cgy * cgz * mg * vel[0] * sin03 / 2.0) + (c5x * l5m * l5x * vel[0] * sin03 / 2.0) - (
                              c5y * l5m * l5z * vel[0] * sin03 / 2.0) + (cgx * l5x * mg * vel[0] * sin03 / 2.0) - (
                              cgy * l5z * mg * vel[0] * sin03 / 2.0) + (c5x * l5m * l3z * vel[0] * cos29 / 2.0) - (
                              c5x * l5m * l4z * vel[1] * cos09 / 2.0) + (c5x * l5m * l4z * vel[2] * cos09 / 2.0) - (
                              c5x * l5m * l4z * vel[3] * cos09 / 2.0) + (c5x * l5m * l4z * vel[4] * cos09 / 2.0) + (
                              cgx * l3z * mg * vel[0] * cos29 / 2.0) - (cgx * l4z * mg * vel[1] * cos09 / 2.0) + (
                              cgx * l4z * mg * vel[2] * cos09 / 2.0) - (cgx * l4z * mg * vel[3] * cos09 / 2.0) + (
                              cgx * l4z * mg * vel[4] * cos09 / 2.0) + (c5y * l5m * l3z * vel[0] * sin32 / 2.0) - (
                              c5y * l5m * l4z * vel[1] * sin15 / 2.0) + (c5y * l5m * l4z * vel[2] * sin15 / 2.0) - (
                              c5y * l5m * l4z * vel[3] * sin15 / 2.0) + (c5y * l5m * l4z * vel[4] * sin15 / 2.0) + (
                              cgy * l3z * mg * vel[0] * sin32 / 2.0) - (cgy * l4z * mg * vel[1] * sin15 / 2.0) + (
                              cgy * l4z * mg * vel[2] * sin15 / 2.0) - (cgy * l4z * mg * vel[3] * sin15 / 2.0) + (
                              cgy * l4z * mg * vel[4] * sin15 / 2.0) - (c5x * c5y * l5m * vel[1] * cos11) + (
                              c5x * c5y * l5m * vel[2] * cos11) - (c5x * c5y * l5m * vel[3] * cos11) - (
                              cgx * cgy * mg * vel[1] * cos11) + (cgx * cgy * mg * vel[2] * cos11) - (
                              cgx * cgy * mg * vel[3] * cos11) - (c5x * c5y * l5m * vel[0] * cos08 / 2.0) + (
                              c5x * c5y * l5m * vel[0] * cos10 / 2.0) + (c5x * c5z * l5m * vel[0] * cos00 / 2.0) - (
                              cgx * cgy * mg * vel[0] * cos08 / 2.0) + (cgx * cgy * mg * vel[0] * cos10 / 2.0) + (
                              cgx * cgz * mg * vel[0] * cos00 / 2.0) - (c5y * l5m * l5x * vel[0] * cos00 / 2.0) + (
                              c5x * l5m * l5z * vel[0] * cos00 / 2.0) - (cgy * l5x * mg * vel[0] * cos00 / 2.0) + (
                              cgx * l5z * mg * vel[0] * cos00 / 2.0) + (c5y * c5z * l5m * vel[0] * sin04 / 2.0) + (
                              cgy * cgz * mg * vel[0] * sin04 / 2.0) + (c5x * l5m * l5x * vel[0] * sin04 / 2.0) + (
                              c5y * l5m * l5z * vel[0] * sin04 / 2.0) + (cgx * l5x * mg * vel[0] * sin04 / 2.0) + (
                              cgy * l5z * mg * vel[0] * sin04 / 2.0) + (c5x * l5m * l4z * vel[0] * cos16 / 2.0) - (
                              c5x * l5m * l3z * vel[1] * cos14 / 2.0) + (c5x * l5m * l3z * vel[2] * cos14 / 2.0) - (
                              c5x * l5m * l3z * vel[3] * cos14 / 2.0) + (c5x * l5m * l3z * vel[4] * cos14 / 2.0) + (
                              cgx * l4z * mg * vel[0] * cos16 / 2.0) - (cgx * l3z * mg * vel[1] * cos14 / 2.0) + (
                              cgx * l3z * mg * vel[2] * cos14 / 2.0) - (cgx * l3z * mg * vel[3] * cos14 / 2.0) + (
                              cgx * l3z * mg * vel[4] * cos14 / 2.0) - (c5y * l5m * l4z * vel[0] * sin21 / 2.0) + (
                              c5y * l5m * l3z * vel[1] * sin16 / 2.0) - (c5y * l5m * l3z * vel[2] * sin16 / 2.0) + (
                              c5y * l5m * l3z * vel[3] * sin16 / 2.0) - (c5y * l5m * l3z * vel[4] * sin16 / 2.0) - (
                              cgy * l4z * mg * vel[0] * sin21 / 2.0) + (cgy * l3z * mg * vel[1] * sin16 / 2.0) - (
                              cgy * l3z * mg * vel[2] * sin16 / 2.0) + (cgy * l3z * mg * vel[3] * sin16 / 2.0) - (
                              cgy * l3z * mg * vel[4] * sin16 / 2.0) + (c5x * l5m * l3z * vel[0] * cos28 / 2.0) + (
                              c5x * l5m * l4z * vel[1] * cos07 / 2.0) - (c5x * l5m * l4z * vel[2] * cos07 / 2.0) + (
                              c5x * l5m * l4z * vel[3] * cos07 / 2.0) + (c5x * l5m * l4z * vel[4] * cos07 / 2.0) + (
                              cgx * l3z * mg * vel[0] * cos28 / 2.0) + (cgx * l4z * mg * vel[1] * cos07 / 2.0) - (
                              cgx * l4z * mg * vel[2] * cos07 / 2.0) + (cgx * l4z * mg * vel[3] * cos07 / 2.0) + (
                              cgx * l4z * mg * vel[4] * cos07 / 2.0) - (c5y * l5m * l3z * vel[0] * sin31 / 2.0) - (
                              c5y * l5m * l4z * vel[1] * sin14 / 2.0) + (c5y * l5m * l4z * vel[2] * sin14 / 2.0) - (
                              c5y * l5m * l4z * vel[3] * sin14 / 2.0) - (c5y * l5m * l4z * vel[4] * sin14 / 2.0) - (
                              cgy * l3z * mg * vel[0] * sin31 / 2.0) - (cgy * l4z * mg * vel[1] * sin14 / 2.0) + (
                              cgy * l4z * mg * vel[2] * sin14 / 2.0) - (cgy * l4z * mg * vel[3] * sin14 / 2.0) - (
                              cgy * l4z * mg * vel[4] * sin14 / 2.0) + (c5x * c5z * l5m * vel[4] * cos02) + (
                              cgx * cgz * mg * vel[4] * cos02) - (c5y * l5m * l5x * vel[1] * cos02) + (
                              c5y * l5m * l5x * vel[2] * cos02) - (c5y * l5m * l5x * vel[3] * cos02) + (
                              c5x * l5m * l5z * vel[4] * cos02) - (cgy * l5x * mg * vel[1] * cos02) + (
                              cgy * l5x * mg * vel[2] * cos02) - (cgy * l5x * mg * vel[3] * cos02) + (
                              cgx * l5z * mg * vel[4] * cos02) + (c5x * l5m * l4z * vel[0] * cos17 / 2.0) + (
                              c5x * l5m * l3z * vel[1] * cos13 / 2.0) - (c5x * l5m * l3z * vel[2] * cos13 / 2.0) + (
                              c5x * l5m * l3z * vel[3] * cos13 / 2.0) + (c5x * l5m * l3z * vel[4] * cos13 / 2.0) - (
                              c5y * c5z * l5m * vel[4] * sin07) + (cgx * l4z * mg * vel[0] * cos17 / 2.0) + (
                              cgx * l3z * mg * vel[1] * cos13 / 2.0) - (cgx * l3z * mg * vel[2] * cos13 / 2.0) + (
                              cgx * l3z * mg * vel[3] * cos13 / 2.0) + (cgx * l3z * mg * vel[4] * cos13 / 2.0) - (
                              cgy * cgz * mg * vel[4] * sin07) - (c5x * l5m * l5x * vel[1] * sin07) + (
                              c5x * l5m * l5x * vel[2] * sin07) - (c5x * l5m * l5x * vel[3] * sin07) - (
                              c5y * l5m * l5z * vel[4] * sin07) - (cgx * l5x * mg * vel[1] * sin07) + (
                              cgx * l5x * mg * vel[2] * sin07) - (cgx * l5x * mg * vel[3] * sin07) - (
                              cgy * l5z * mg * vel[4] * sin07) + (c5x * c5z * l5m * vel[0] * cos01 / 2.0) + (
                              cgx * cgz * mg * vel[0] * cos01 / 2.0) + (c5y * l5m * l4z * vel[0] * sin22 / 2.0) + (
                              c5y * l5m * l3z * vel[1] * sin17 / 2.0) - (c5y * l5m * l3z * vel[2] * sin17 / 2.0) + (
                              c5y * l5m * l3z * vel[3] * sin17 / 2.0) + (c5y * l5m * l3z * vel[4] * sin17 / 2.0) + (
                              cgy * l4z * mg * vel[0] * sin22 / 2.0) + (cgy * l3z * mg * vel[1] * sin17 / 2.0) - (
                              cgy * l3z * mg * vel[2] * sin17 / 2.0) + (cgy * l3z * mg * vel[3] * sin17 / 2.0) + (
                              cgy * l3z * mg * vel[4] * sin17 / 2.0)
        C[2, 0] = (vel[0] * ((I3x * sin18) - (I3z * sin18) + (I5x * ((sin10 / 2.0) + (sin09 / 2.0) + sin05) / 2.0) + (
                    Igx * ((sin10 / 2.0) + (sin09 / 2.0) + sin05) / 2.0) - (
                                         I5y * ((sin10 / 2.0) + (sin09 / 2.0) - sin05) / 2.0) - (
                                         Igy * ((sin10 / 2.0) + (sin09 / 2.0) - sin05) / 2.0) + (I4x * sin05) - (
                                         I4z * sin05) - (I5z * sin05) - (Igz * sin05) - (l5m * (
                    (c5z * sin04) - (l5x * cos00) + (l5z * sin04) + (2.0 * c5y * sin06) - (l4z * sin21) + (
                        l4z * sin22) + (l5x * cos01) - (c5z * sin03) - (l5z * sin03)) * ((c5z * cos00) + (
                    l5z * cos00) + (l5x * sin04) + (2.0 * c5y * cos04) - (l4z * cos16) - (l3z * cos28) + (
                                                                                                     l4z * cos17) - (
                                                                                                     2.0 * l2x * sin07) - (
                                                                                                     c5z * cos01) - (
                                                                                                     l5z * cos01) - (
                                                                                                     l5x * sin03) + (
                                                                                                     l3z * cos29)) / 2.0) - (
                                         mg * ((cgz * sin04) - (l5x * cos00) + (l5z * sin04) + (2.0 * cgy * sin06) - (
                                             l4z * sin21) + (l4z * sin22) + (l5x * cos01) - (cgz * sin03) - (
                                                           l5z * sin03)) * (
                                                     (cgz * cos00) + (l5z * cos00) + (l5x * sin04) + (
                                                         2.0 * cgy * cos04) - (l4z * cos16) - (l3z * cos28) + (
                                                                 l4z * cos17) - (2.0 * l2x * sin07) - (cgz * cos01) - (
                                                                 l5z * cos01) - (l5x * sin03) + (
                                                                 l3z * cos29)) / 2.0) + (l5m * (
                    (c5z * cos00) + (l5z * cos00) + (l5x * sin04) + (l4z * cos16) + (2.0 * c5x * sin06) + (
                        l4z * cos17) + (c5z * cos01) + (l5z * cos01) + (l5x * sin03)) * ((c5z * sin04) - (
                    l5x * cos00) + (l5z * sin04) - (2.0 * c5x * cos04) + (l4z * sin21) + (l3z * sin31) - (
                                                                                                     2.0 * l2x * cos02) + (
                                                                                                     l4z * sin22) - (
                                                                                                     l5x * cos01) + (
                                                                                                     c5z * sin03) + (
                                                                                                     l5z * sin03) + (
                                                                                                     l3z * sin32)) / 2.0) + (
                                         mg * ((cgz * cos00) + (l5z * cos00) + (l5x * sin04) + (l4z * cos16) + (
                                             2.0 * cgx * sin06) + (l4z * cos17) + (cgz * cos01) + (l5z * cos01) + (
                                                           l5x * sin03)) * (
                                                     (cgz * sin04) - (l5x * cos00) + (l5z * sin04) - (
                                                         2.0 * cgx * cos04) + (l4z * sin21) + (l3z * sin31) - (
                                                                 2.0 * l2x * cos02) + (l4z * sin22) - (l5x * cos01) + (
                                                                 cgz * sin03) + (l5z * sin03) + (
                                                                 l3z * sin32)) / 2.0) + (
                                         2.0 * l4m * ((c4z * cos04) + (c4x * sin06) + (l4z * cos18)) * (
                                             (c4z * sin06) - (c4x * cos04) - l2x + (l3z * sin33) + (l4z * sin30))) + (
                                         l5m * ((c5x * cos00) + (c5y * sin04) - (c5x * cos01) + (c5y * sin03)) * (
                                             (c5y * cos00) - (c5x * sin04) + (c5y * cos01) + (c5x * sin03)) / 2.0) + (
                                         mg * ((cgx * cos00) + (cgy * sin04) - (cgx * cos01) + (cgy * sin03)) * (
                                             (cgy * cos00) - (cgx * sin04) + (cgy * cos01) + (cgx * sin03)) / 2.0) - (
                                         2.0 * l3m * ((c3z * cos18) + (c3x * sin30)) * (
                                             l2x - (l3z * sin33) + (c3x * cos18) - (c3z * sin30)))) / 2.0) - (vel[4] * (
                    (I5y * sin01 / 2.0) - (I5x * sin02 / 2.0) - (I5x * sin01 / 2.0) + (I5y * sin02 / 2.0) - (
                        Igx * sin01 / 2.0) - (Igx * sin02 / 2.0) + (Igy * sin01 / 2.0) + (Igy * sin02 / 2.0) + (
                                I5z * sin06) + (Igz * sin06) + (pow(c5x, 2.0) * l5m * sin06) + (
                                pow(c5y, 2.0) * l5m * sin06) + (pow(cgx, 2.0) * mg * sin06) + (
                                pow(cgy, 2.0) * mg * sin06) + (pow(c5x, 2.0) * l5m * sin01 / 2.0) + (
                                pow(c5x, 2.0) * l5m * sin02 / 2.0) - (pow(c5y, 2.0) * l5m * sin01 / 2.0) - (
                                pow(c5y, 2.0) * l5m * sin02 / 2.0) + (pow(cgx, 2.0) * mg * sin01 / 2.0) + (
                                pow(cgx, 2.0) * mg * sin02 / 2.0) - (pow(cgy, 2.0) * mg * sin01 / 2.0) - (
                                pow(cgy, 2.0) * mg * sin02 / 2.0) + (c5x * l5m * l4z * cos16) + (
                                cgx * l4z * mg * cos16) - (c5y * l5m * l4z * sin21) - (cgy * l4z * mg * sin21) + (
                                c5x * l5m * l4z * cos17) + (cgx * l4z * mg * cos17) + (c5x * c5z * l5m * cos01) + (
                                cgx * cgz * mg * cos01) + (c5y * l5m * l4z * sin22) + (cgy * l4z * mg * sin22) + (
                                c5y * l5m * l5x * cos01) + (c5x * l5m * l5z * cos01) + (cgy * l5x * mg * cos01) + (
                                cgx * l5z * mg * cos01) - (c5y * c5z * l5m * sin03) - (cgy * cgz * mg * sin03) + (
                                c5x * l5m * l5x * sin03) - (c5y * l5m * l5z * sin03) + (cgx * l5x * mg * sin03) - (
                                cgy * l5z * mg * sin03) - (c5x * c5y * l5m * cos08) + (c5x * c5y * l5m * cos10) + (
                                c5x * c5z * l5m * cos00) - (cgx * cgy * mg * cos08) + (cgx * cgy * mg * cos10) + (
                                cgx * cgz * mg * cos00) - (c5y * l5m * l5x * cos00) + (c5x * l5m * l5z * cos00) - (
                                cgy * l5x * mg * cos00) + (cgx * l5z * mg * cos00) + (c5y * c5z * l5m * sin04) + (
                                cgy * cgz * mg * sin04) + (c5x * l5m * l5x * sin04) + (c5y * l5m * l5z * sin04) + (
                                cgx * l5x * mg * sin04) + (cgy * l5z * mg * sin04)) / 2.0)
        C[2, 1] = (I5y * vel[4] * sin00 / 2.0) - (I5x * vel[4] * sin00 / 2.0) - (Igx * vel[4] * sin00 / 2.0) + (
                    Igy * vel[4] * sin00 / 2.0) + (pow(c5x, 2.0) * l5m * vel[4] * sin00 / 2.0) - (
                              pow(c5y, 2.0) * l5m * vel[4] * sin00 / 2.0) + (
                              pow(cgx, 2.0) * mg * vel[4] * sin00 / 2.0) - (
                              pow(cgy, 2.0) * mg * vel[4] * sin00 / 2.0) + (l3z * l4z * mg * vel[1] * sin24) + (
                              l4z * l5z * mg * vel[3] * sin08) + (c4x * l4m * l3z * vel[1] * cos12) - (
                              c5x * l5m * l4z * vel[3] * cos09 / 2.0) + (c5x * l5m * l4z * vel[4] * cos09 / 2.0) - (
                              cgx * l4z * mg * vel[3] * cos09 / 2.0) + (cgx * l4z * mg * vel[4] * cos09 / 2.0) + (
                              l5m * l5x * l3z * vel[1] * cos12) + (l5x * l3z * mg * vel[1] * cos12) - (
                              c5y * l5m * l4z * vel[3] * sin15 / 2.0) + (c5y * l5m * l4z * vel[4] * sin15 / 2.0) + (
                              c4z * l4m * l3z * vel[1] * sin11) + (c5z * l5m * l3z * vel[1] * sin11) - (
                              cgy * l4z * mg * vel[3] * sin15 / 2.0) + (cgy * l4z * mg * vel[4] * sin15 / 2.0) + (
                              cgz * l3z * mg * vel[1] * sin11) + (l5m * l3z * l5z * vel[1] * sin11) + (
                              l3z * l5z * mg * vel[1] * sin11) + (c5x * c5y * l5m * vel[4] * cos11) + (
                              cgx * cgy * mg * vel[4] * cos11) + (c5x * l5m * l3z * vel[1] * cos14 / 2.0) + (
                              cgx * l3z * mg * vel[1] * cos14 / 2.0) - (c5y * l5m * l3z * vel[1] * sin16 / 2.0) - (
                              cgy * l3z * mg * vel[1] * sin16 / 2.0) - (c5x * l5m * l4z * vel[3] * cos07 / 2.0) - (
                              c5x * l5m * l4z * vel[4] * cos07 / 2.0) - (cgx * l4z * mg * vel[3] * cos07 / 2.0) - (
                              cgx * l4z * mg * vel[4] * cos07 / 2.0) + (c5y * l5m * l4z * vel[3] * sin14 / 2.0) + (
                              c5y * l5m * l4z * vel[4] * sin14 / 2.0) + (cgy * l4z * mg * vel[3] * sin14 / 2.0) + (
                              cgy * l4z * mg * vel[4] * sin14 / 2.0) + (c5y * l5m * l5x * vel[4] * cos02) + (
                              c3x * l3m * l3z * vel[1] * cos30) - (c4x * l4m * l4z * vel[3] * cos03) + (
                              cgy * l5x * mg * vel[4] * cos02) + (c5x * l5m * l3z * vel[1] * cos13 / 2.0) + (
                              cgx * l3z * mg * vel[1] * cos13 / 2.0) - (l5m * l5x * l4z * vel[3] * cos03) - (
                              l5x * l4z * mg * vel[3] * cos03) + (c5x * l5m * l5x * vel[4] * sin07) + (
                              c3z * l3m * l3z * vel[1] * sin24) + (c4z * l4m * l4z * vel[3] * sin08) + (
                              c5z * l5m * l4z * vel[3] * sin08) + (cgx * l5x * mg * vel[4] * sin07) + (
                              cgz * l4z * mg * vel[3] * sin08) + (c5y * l5m * l3z * vel[1] * sin17 / 2.0) + (
                              cgy * l3z * mg * vel[1] * sin17 / 2.0) + (l4m * l3z * l4z * vel[1] * sin24) + (
                              l5m * l3z * l4z * vel[1] * sin24) + (l5m * l4z * l5z * vel[3] * sin08)
        C[2, 2] = -(vel[4] * ((I5y * sin00) - (I5x * sin00) - (Igx * sin00) + (Igy * sin00) + (
                    2.0 * l5m * ((c5y * cos02) + (c5x * sin07)) * (
                        l5x + (c5x * cos02) - (c5y * sin07) + (l4z * sin08))) + (
                                          2.0 * mg * ((cgy * cos02) + (cgx * sin07)) * (
                                              l5x + (cgx * cos02) - (cgy * sin07) + (l4z * sin08)))) / 2.0) - (
                              l4z * vel[3] * ((c4z * l4m * sin08) - (l5m * l5x * cos03) - (l5x * mg * cos03) - (
                                  c4x * l4m * cos03) + (c5z * l5m * sin08) + (cgz * mg * sin08) + (
                                                          l5m * l5z * sin08) + (l5z * mg * sin08) - (
                                                          c5x * l5m * cos03 * cos02) - (cgx * mg * cos03 * cos02) + (
                                                          c5y * l5m * cos03 * sin07) + (cgy * mg * cos03 * sin07)))
        C[2, 3] = (I5y * vel[4] * sin00 / 2.0) - (I5x * vel[4] * sin00 / 2.0) - (Igx * vel[4] * sin00 / 2.0) + (
                    Igy * vel[4] * sin00 / 2.0) - (c5x * c5y * l5m * vel[4]) - (cgx * cgy * mg * vel[4]) + (
                              pow(c5x, 2.0) * l5m * vel[4] * sin00 / 2.0) - (
                              pow(c5y, 2.0) * l5m * vel[4] * sin00 / 2.0) + (
                              pow(cgx, 2.0) * mg * vel[4] * sin00 / 2.0) - (
                              pow(cgy, 2.0) * mg * vel[4] * sin00 / 2.0) + (l4z * l5z * mg * vel[1] * sin08) - (
                              l4z * l5z * mg * vel[2] * sin08) + (l4z * l5z * mg * vel[3] * sin08) + (
                              2.0 * c5x * c5y * l5m * vel[4] * pow(cos02, 2.0)) + (
                              2.0 * cgx * cgy * mg * vel[4] * pow(cos02, 2.0)) + (c5y * l5m * l5x * vel[4] * cos02) - (
                              c4x * l4m * l4z * vel[1] * cos03) + (c4x * l4m * l4z * vel[2] * cos03) - (
                              c4x * l4m * l4z * vel[3] * cos03) + (cgy * l5x * mg * vel[4] * cos02) - (
                              l5m * l5x * l4z * vel[1] * cos03) + (l5m * l5x * l4z * vel[2] * cos03) - (
                              l5m * l5x * l4z * vel[3] * cos03) - (l5x * l4z * mg * vel[1] * cos03) + (
                              l5x * l4z * mg * vel[2] * cos03) - (l5x * l4z * mg * vel[3] * cos03) + (
                              c5x * l5m * l5x * vel[4] * sin07) + (c4z * l4m * l4z * vel[1] * sin08) - (
                              c4z * l4m * l4z * vel[2] * sin08) + (c4z * l4m * l4z * vel[3] * sin08) + (
                              c5z * l5m * l4z * vel[1] * sin08) - (c5z * l5m * l4z * vel[2] * sin08) + (
                              c5z * l5m * l4z * vel[3] * sin08) + (cgx * l5x * mg * vel[4] * sin07) + (
                              cgz * l4z * mg * vel[1] * sin08) - (cgz * l4z * mg * vel[2] * sin08) + (
                              cgz * l4z * mg * vel[3] * sin08) + (l5m * l4z * l5z * vel[1] * sin08) - (
                              l5m * l4z * l5z * vel[2] * sin08) + (l5m * l4z * l5z * vel[3] * sin08) - (
                              c5x * l5m * l4z * vel[1] * cos03 * cos02) + (c5x * l5m * l4z * vel[2] * cos03 * cos02) - (
                              c5x * l5m * l4z * vel[3] * cos03 * cos02) - (cgx * l4z * mg * vel[1] * cos03 * cos02) + (
                              cgx * l4z * mg * vel[2] * cos03 * cos02) - (cgx * l4z * mg * vel[3] * cos03 * cos02) + (
                              c5y * l5m * l4z * vel[1] * cos03 * sin07) - (c5y * l5m * l4z * vel[2] * cos03 * sin07) + (
                              c5y * l5m * l4z * vel[3] * cos03 * sin07) + (c5y * l5m * l4z * vel[4] * cos02 * sin08) + (
                              cgy * l4z * mg * vel[1] * cos03 * sin07) - (cgy * l4z * mg * vel[2] * cos03 * sin07) + (
                              cgy * l4z * mg * vel[3] * cos03 * sin07) + (cgy * l4z * mg * vel[4] * cos02 * sin08) + (
                              c5x * l5m * l4z * vel[4] * sin08 * sin07) + (cgx * l4z * mg * vel[4] * sin08 * sin07)
        C[2, 4] = (I5x * vel[2] * sin00 / 2.0) - (I5x * vel[1] * sin00 / 2.0) - (I5x * vel[3] * sin00 / 2.0) + (
                    I5y * vel[1] * sin00 / 2.0) - (I5y * vel[2] * sin00 / 2.0) + (I5y * vel[3] * sin00 / 2.0) - (
                              Igx * vel[1] * sin00 / 2.0) + (Igx * vel[2] * sin00 / 2.0) - (
                              Igx * vel[3] * sin00 / 2.0) + (Igy * vel[1] * sin00 / 2.0) - (
                              Igy * vel[2] * sin00 / 2.0) + (Igy * vel[3] * sin00 / 2.0) + (
                              I5x * vel[0] * sin01 / 4.0) + (I5x * vel[0] * sin02 / 4.0) - (
                              I5y * vel[0] * sin01 / 4.0) - (I5y * vel[0] * sin02 / 4.0) + (
                              Igx * vel[0] * sin01 / 4.0) + (Igx * vel[0] * sin02 / 4.0) - (
                              Igy * vel[0] * sin01 / 4.0) - (Igy * vel[0] * sin02 / 4.0) - (
                              I5z * vel[0] * sin06 / 2.0) - (Igz * vel[0] * sin06 / 2.0) + (
                              pow(c5x, 2.0) * l5m * vel[1] * sin00 / 2.0) - (
                              pow(c5x, 2.0) * l5m * vel[2] * sin00 / 2.0) + (
                              pow(c5x, 2.0) * l5m * vel[3] * sin00 / 2.0) - (
                              pow(c5y, 2.0) * l5m * vel[1] * sin00 / 2.0) + (
                              pow(c5y, 2.0) * l5m * vel[2] * sin00 / 2.0) - (
                              pow(c5y, 2.0) * l5m * vel[3] * sin00 / 2.0) + (
                              pow(cgx, 2.0) * mg * vel[1] * sin00 / 2.0) - (
                              pow(cgx, 2.0) * mg * vel[2] * sin00 / 2.0) + (
                              pow(cgx, 2.0) * mg * vel[3] * sin00 / 2.0) - (
                              pow(cgy, 2.0) * mg * vel[1] * sin00 / 2.0) + (
                              pow(cgy, 2.0) * mg * vel[2] * sin00 / 2.0) - (
                              pow(cgy, 2.0) * mg * vel[3] * sin00 / 2.0) - (
                              pow(c5x, 2.0) * l5m * vel[0] * sin01 / 4.0) - (
                              pow(c5x, 2.0) * l5m * vel[0] * sin02 / 4.0) + (
                              pow(c5y, 2.0) * l5m * vel[0] * sin01 / 4.0) + (
                              pow(c5y, 2.0) * l5m * vel[0] * sin02 / 4.0) - (
                              pow(cgx, 2.0) * mg * vel[0] * sin01 / 4.0) - (
                              pow(cgx, 2.0) * mg * vel[0] * sin02 / 4.0) + (
                              pow(cgy, 2.0) * mg * vel[0] * sin01 / 4.0) + (
                              pow(cgy, 2.0) * mg * vel[0] * sin02 / 4.0) - (
                              pow(c5x, 2.0) * l5m * vel[0] * sin06 / 2.0) - (
                              pow(c5y, 2.0) * l5m * vel[0] * sin06 / 2.0) - (
                              pow(cgx, 2.0) * mg * vel[0] * sin06 / 2.0) - (
                              pow(cgy, 2.0) * mg * vel[0] * sin06 / 2.0) - (c5y * l5m * l5x * vel[0] * cos01 / 2.0) - (
                              c5x * l5m * l5z * vel[0] * cos01 / 2.0) - (cgy * l5x * mg * vel[0] * cos01 / 2.0) - (
                              cgx * l5z * mg * vel[0] * cos01 / 2.0) + (c5y * c5z * l5m * vel[0] * sin03 / 2.0) + (
                              cgy * cgz * mg * vel[0] * sin03 / 2.0) - (c5x * l5m * l5x * vel[0] * sin03 / 2.0) + (
                              c5y * l5m * l5z * vel[0] * sin03 / 2.0) - (cgx * l5x * mg * vel[0] * sin03 / 2.0) + (
                              cgy * l5z * mg * vel[0] * sin03 / 2.0) + (c5x * l5m * l4z * vel[1] * cos09 / 2.0) - (
                              c5x * l5m * l4z * vel[2] * cos09 / 2.0) + (c5x * l5m * l4z * vel[3] * cos09 / 2.0) - (
                              c5x * l5m * l4z * vel[4] * cos09 / 2.0) + (cgx * l4z * mg * vel[1] * cos09 / 2.0) - (
                              cgx * l4z * mg * vel[2] * cos09 / 2.0) + (cgx * l4z * mg * vel[3] * cos09 / 2.0) - (
                              cgx * l4z * mg * vel[4] * cos09 / 2.0) + (c5y * l5m * l4z * vel[1] * sin15 / 2.0) - (
                              c5y * l5m * l4z * vel[2] * sin15 / 2.0) + (c5y * l5m * l4z * vel[3] * sin15 / 2.0) - (
                              c5y * l5m * l4z * vel[4] * sin15 / 2.0) + (cgy * l4z * mg * vel[1] * sin15 / 2.0) - (
                              cgy * l4z * mg * vel[2] * sin15 / 2.0) + (cgy * l4z * mg * vel[3] * sin15 / 2.0) - (
                              cgy * l4z * mg * vel[4] * sin15 / 2.0) + (c5x * c5y * l5m * vel[1] * cos11) - (
                              c5x * c5y * l5m * vel[2] * cos11) + (c5x * c5y * l5m * vel[3] * cos11) + (
                              cgx * cgy * mg * vel[1] * cos11) - (cgx * cgy * mg * vel[2] * cos11) + (
                              cgx * cgy * mg * vel[3] * cos11) + (c5x * c5y * l5m * vel[0] * cos08 / 2.0) - (
                              c5x * c5y * l5m * vel[0] * cos10 / 2.0) - (c5x * c5z * l5m * vel[0] * cos00 / 2.0) + (
                              cgx * cgy * mg * vel[0] * cos08 / 2.0) - (cgx * cgy * mg * vel[0] * cos10 / 2.0) - (
                              cgx * cgz * mg * vel[0] * cos00 / 2.0) + (c5y * l5m * l5x * vel[0] * cos00 / 2.0) - (
                              c5x * l5m * l5z * vel[0] * cos00 / 2.0) + (cgy * l5x * mg * vel[0] * cos00 / 2.0) - (
                              cgx * l5z * mg * vel[0] * cos00 / 2.0) - (c5y * c5z * l5m * vel[0] * sin04 / 2.0) - (
                              cgy * cgz * mg * vel[0] * sin04 / 2.0) - (c5x * l5m * l5x * vel[0] * sin04 / 2.0) - (
                              c5y * l5m * l5z * vel[0] * sin04 / 2.0) - (cgx * l5x * mg * vel[0] * sin04 / 2.0) - (
                              cgy * l5z * mg * vel[0] * sin04 / 2.0) - (c5x * l5m * l4z * vel[0] * cos16 / 2.0) - (
                              cgx * l4z * mg * vel[0] * cos16 / 2.0) + (c5y * l5m * l4z * vel[0] * sin21 / 2.0) + (
                              cgy * l4z * mg * vel[0] * sin21 / 2.0) - (c5x * l5m * l4z * vel[1] * cos07 / 2.0) + (
                              c5x * l5m * l4z * vel[2] * cos07 / 2.0) - (c5x * l5m * l4z * vel[3] * cos07 / 2.0) - (
                              c5x * l5m * l4z * vel[4] * cos07 / 2.0) - (cgx * l4z * mg * vel[1] * cos07 / 2.0) + (
                              cgx * l4z * mg * vel[2] * cos07 / 2.0) - (cgx * l4z * mg * vel[3] * cos07 / 2.0) - (
                              cgx * l4z * mg * vel[4] * cos07 / 2.0) + (c5y * l5m * l4z * vel[1] * sin14 / 2.0) - (
                              c5y * l5m * l4z * vel[2] * sin14 / 2.0) + (c5y * l5m * l4z * vel[3] * sin14 / 2.0) + (
                              c5y * l5m * l4z * vel[4] * sin14 / 2.0) + (cgy * l4z * mg * vel[1] * sin14 / 2.0) - (
                              cgy * l4z * mg * vel[2] * sin14 / 2.0) + (cgy * l4z * mg * vel[3] * sin14 / 2.0) + (
                              cgy * l4z * mg * vel[4] * sin14 / 2.0) - (c5x * c5z * l5m * vel[4] * cos02) - (
                              cgx * cgz * mg * vel[4] * cos02) + (c5y * l5m * l5x * vel[1] * cos02) - (
                              c5y * l5m * l5x * vel[2] * cos02) + (c5y * l5m * l5x * vel[3] * cos02) - (
                              c5x * l5m * l5z * vel[4] * cos02) + (cgy * l5x * mg * vel[1] * cos02) - (
                              cgy * l5x * mg * vel[2] * cos02) + (cgy * l5x * mg * vel[3] * cos02) - (
                              cgx * l5z * mg * vel[4] * cos02) - (c5x * l5m * l4z * vel[0] * cos17 / 2.0) + (
                              c5y * c5z * l5m * vel[4] * sin07) - (cgx * l4z * mg * vel[0] * cos17 / 2.0) + (
                              cgy * cgz * mg * vel[4] * sin07) + (c5x * l5m * l5x * vel[1] * sin07) - (
                              c5x * l5m * l5x * vel[2] * sin07) + (c5x * l5m * l5x * vel[3] * sin07) + (
                              c5y * l5m * l5z * vel[4] * sin07) + (cgx * l5x * mg * vel[1] * sin07) - (
                              cgx * l5x * mg * vel[2] * sin07) + (cgx * l5x * mg * vel[3] * sin07) + (
                              cgy * l5z * mg * vel[4] * sin07) - (c5x * c5z * l5m * vel[0] * cos01 / 2.0) - (
                              cgx * cgz * mg * vel[0] * cos01 / 2.0) - (c5y * l5m * l4z * vel[0] * sin22 / 2.0) - (
                              cgy * l4z * mg * vel[0] * sin22 / 2.0)
        C[3, 0] = (vel[4] * ((I5y * sin01 / 2.0) - (I5x * sin02 / 2.0) - (I5x * sin01 / 2.0) + (I5y * sin02 / 2.0) - (
                    Igx * sin01 / 2.0) - (Igx * sin02 / 2.0) + (Igy * sin01 / 2.0) + (Igy * sin02 / 2.0) + (
                                         I5z * sin06) + (Igz * sin06) + (pow(c5x, 2.0) * l5m * sin06) + (
                                         pow(c5y, 2.0) * l5m * sin06) + (pow(cgx, 2.0) * mg * sin06) + (
                                         pow(cgy, 2.0) * mg * sin06) + (pow(c5x, 2.0) * l5m * sin01 / 2.0) + (
                                         pow(c5x, 2.0) * l5m * sin02 / 2.0) - (pow(c5y, 2.0) * l5m * sin01 / 2.0) - (
                                         pow(c5y, 2.0) * l5m * sin02 / 2.0) + (pow(cgx, 2.0) * mg * sin01 / 2.0) + (
                                         pow(cgx, 2.0) * mg * sin02 / 2.0) - (pow(cgy, 2.0) * mg * sin01 / 2.0) - (
                                         pow(cgy, 2.0) * mg * sin02 / 2.0) + (c5x * c5z * l5m * cos01) + (
                                         cgx * cgz * mg * cos01) + (c5y * l5m * l5x * cos01) + (
                                         c5x * l5m * l5z * cos01) + (cgy * l5x * mg * cos01) + (
                                         cgx * l5z * mg * cos01) - (c5y * c5z * l5m * sin03) - (
                                         cgy * cgz * mg * sin03) + (c5x * l5m * l5x * sin03) - (
                                         c5y * l5m * l5z * sin03) + (cgx * l5x * mg * sin03) - (
                                         cgy * l5z * mg * sin03) - (c5x * c5y * l5m * cos08) + (
                                         c5x * c5y * l5m * cos10) + (c5x * c5z * l5m * cos00) - (
                                         cgx * cgy * mg * cos08) + (cgx * cgy * mg * cos10) + (
                                         cgx * cgz * mg * cos00) - (c5y * l5m * l5x * cos00) + (
                                         c5x * l5m * l5z * cos00) - (cgy * l5x * mg * cos00) + (
                                         cgx * l5z * mg * cos00) + (c5y * c5z * l5m * sin04) + (
                                         cgy * cgz * mg * sin04) + (c5x * l5m * l5x * sin04) + (
                                         c5y * l5m * l5z * sin04) + (cgx * l5x * mg * sin04) + (
                                         cgy * l5z * mg * sin04)) / 2.0) - (vel[0] * (
                    (I5x * (cos00 + cos01) * (sin03 + sin04) / 2.0) + (
                        Igx * (cos00 + cos01) * (sin03 + sin04) / 2.0) + (
                                2.0 * l4m * ((c4z * cos04) + (c4x * sin06)) * (
                                    (c4z * sin06) - (c4x * cos04) - l2x + (l3z * sin33) + (l4z * sin30))) - (l5m * (
                        (c5z * sin04) - (l5x * cos00) + (l5z * sin04) + (2.0 * c5y * sin06) + (l5x * cos01) - (
                            c5z * sin03) - (l5z * sin03)) * ((c5z * cos00) + (l5z * cos00) + (l5x * sin04) + (
                        2.0 * c5y * cos04) - (l4z * cos16) - (l3z * cos28) + (l4z * cos17) - (2.0 * l2x * sin07) - (
                                                                         c5z * cos01) - (l5z * cos01) - (
                                                                         l5x * sin03) + (l3z * cos29)) / 2.0) - (mg * (
                        (cgz * sin04) - (l5x * cos00) + (l5z * sin04) + (2.0 * cgy * sin06) + (l5x * cos01) - (
                            cgz * sin03) - (l5z * sin03)) * ((cgz * cos00) + (l5z * cos00) + (l5x * sin04) + (
                        2.0 * cgy * cos04) - (l4z * cos16) - (l3z * cos28) + (l4z * cos17) - (2.0 * l2x * sin07) - (
                                                                         cgz * cos01) - (l5z * cos01) - (
                                                                         l5x * sin03) + (l3z * cos29)) / 2.0) + (l5m * (
                        (c5z * cos00) + (l5z * cos00) + (l5x * sin04) + (2.0 * c5x * sin06) + (c5z * cos01) + (
                            l5z * cos01) + (l5x * sin03)) * ((c5z * sin04) - (l5x * cos00) + (l5z * sin04) - (
                        2.0 * c5x * cos04) + (l4z * sin21) + (l3z * sin31) - (2.0 * l2x * cos02) + (l4z * sin22) - (
                                                                         l5x * cos01) + (c5z * sin03) + (
                                                                         l5z * sin03) + (l3z * sin32)) / 2.0) + (mg * (
                        (cgz * cos00) + (l5z * cos00) + (l5x * sin04) + (2.0 * cgx * sin06) + (cgz * cos01) + (
                            l5z * cos01) + (l5x * sin03)) * ((cgz * sin04) - (l5x * cos00) + (l5z * sin04) - (
                        2.0 * cgx * cos04) + (l4z * sin21) + (l3z * sin31) - (2.0 * l2x * cos02) + (l4z * sin22) - (
                                                                         l5x * cos01) + (cgz * sin03) + (
                                                                         l5z * sin03) + (l3z * sin32)) / 2.0) + (
                                2.0 * I4x * cos04 * sin06) - (2.0 * I4z * cos04 * sin06) - (
                                2.0 * I5z * cos04 * sin06) - (2.0 * Igz * cos04 * sin06) + (
                                I5y * (cos00 - cos01) * (sin03 - sin04) / 2.0) + (
                                Igy * (cos00 - cos01) * (sin03 - sin04) / 2.0) + (
                                l5m * ((c5x * cos00) + (c5y * sin04) - (c5x * cos01) + (c5y * sin03)) * (
                                    (c5y * cos00) - (c5x * sin04) + (c5y * cos01) + (c5x * sin03)) / 2.0) + (
                                mg * ((cgx * cos00) + (cgy * sin04) - (cgx * cos01) + (cgy * sin03)) * (
                                    (cgy * cos00) - (cgx * sin04) + (cgy * cos01) + (cgx * sin03)) / 2.0)) / 2.0)
        C[3, 1] = (I5x * vel[4] * sin00 / 2.0) - (I5y * vel[4] * sin00 / 2.0) + (Igx * vel[4] * sin00 / 2.0) - (
                    Igy * vel[4] * sin00 / 2.0) - (pow(c5x, 2.0) * l5m * vel[4] * sin00 / 2.0) + (
                              pow(c5y, 2.0) * l5m * vel[4] * sin00 / 2.0) - (
                              pow(cgx, 2.0) * mg * vel[4] * sin00 / 2.0) + (
                              pow(cgy, 2.0) * mg * vel[4] * sin00 / 2.0) + (l4z * l5z * mg * vel[1] * sin08) - (
                              l4z * l5z * mg * vel[2] * sin08) - (c4x * l4m * l3z * vel[1] * cos12) - (
                              c5x * l5m * l4z * vel[1] * cos09 / 2.0) + (c5x * l5m * l4z * vel[2] * cos09 / 2.0) - (
                              cgx * l4z * mg * vel[1] * cos09 / 2.0) + (cgx * l4z * mg * vel[2] * cos09 / 2.0) - (
                              l5m * l5x * l3z * vel[1] * cos12) - (l5x * l3z * mg * vel[1] * cos12) - (
                              c5y * l5m * l4z * vel[1] * sin15 / 2.0) + (c5y * l5m * l4z * vel[2] * sin15 / 2.0) - (
                              c4z * l4m * l3z * vel[1] * sin11) - (c5z * l5m * l3z * vel[1] * sin11) - (
                              cgy * l4z * mg * vel[1] * sin15 / 2.0) + (cgy * l4z * mg * vel[2] * sin15 / 2.0) - (
                              cgz * l3z * mg * vel[1] * sin11) - (l5m * l3z * l5z * vel[1] * sin11) - (
                              l3z * l5z * mg * vel[1] * sin11) - (c5x * c5y * l5m * vel[4] * cos11) - (
                              cgx * cgy * mg * vel[4] * cos11) - (c5x * l5m * l3z * vel[1] * cos14 / 2.0) - (
                              cgx * l3z * mg * vel[1] * cos14 / 2.0) + (c5y * l5m * l3z * vel[1] * sin16 / 2.0) + (
                              cgy * l3z * mg * vel[1] * sin16 / 2.0) - (c5x * l5m * l4z * vel[1] * cos07 / 2.0) + (
                              c5x * l5m * l4z * vel[2] * cos07 / 2.0) - (cgx * l4z * mg * vel[1] * cos07 / 2.0) + (
                              cgx * l4z * mg * vel[2] * cos07 / 2.0) + (c5y * l5m * l4z * vel[1] * sin14 / 2.0) - (
                              c5y * l5m * l4z * vel[2] * sin14 / 2.0) + (cgy * l4z * mg * vel[1] * sin14 / 2.0) - (
                              cgy * l4z * mg * vel[2] * sin14 / 2.0) - (c5y * l5m * l5x * vel[4] * cos02) - (
                              c4x * l4m * l4z * vel[1] * cos03) + (c4x * l4m * l4z * vel[2] * cos03) - (
                              cgy * l5x * mg * vel[4] * cos02) - (c5x * l5m * l3z * vel[1] * cos13 / 2.0) - (
                              cgx * l3z * mg * vel[1] * cos13 / 2.0) - (l5m * l5x * l4z * vel[1] * cos03) + (
                              l5m * l5x * l4z * vel[2] * cos03) - (l5x * l4z * mg * vel[1] * cos03) + (
                              l5x * l4z * mg * vel[2] * cos03) - (c5x * l5m * l5x * vel[4] * sin07) + (
                              c4z * l4m * l4z * vel[1] * sin08) - (c4z * l4m * l4z * vel[2] * sin08) + (
                              c5z * l5m * l4z * vel[1] * sin08) - (c5z * l5m * l4z * vel[2] * sin08) - (
                              cgx * l5x * mg * vel[4] * sin07) + (cgz * l4z * mg * vel[1] * sin08) - (
                              cgz * l4z * mg * vel[2] * sin08) - (c5y * l5m * l3z * vel[1] * sin17 / 2.0) - (
                              cgy * l3z * mg * vel[1] * sin17 / 2.0) + (l5m * l4z * l5z * vel[1] * sin08) - (
                              l5m * l4z * l5z * vel[2] * sin08)
        C[3, 2] = (I5y * vel[4] * sin00 / 2.0) - (I5x * vel[4] * sin00 / 2.0) - (Igx * vel[4] * sin00 / 2.0) + (
                    Igy * vel[4] * sin00 / 2.0) - (c5x * c5y * l5m * vel[4]) - (cgx * cgy * mg * vel[4]) + (
                              pow(c5x, 2.0) * l5m * vel[4] * sin00 / 2.0) - (
                              pow(c5y, 2.0) * l5m * vel[4] * sin00 / 2.0) + (
                              pow(cgx, 2.0) * mg * vel[4] * sin00 / 2.0) - (
                              pow(cgy, 2.0) * mg * vel[4] * sin00 / 2.0) - (l4z * l5z * mg * vel[1] * sin08) + (
                              l4z * l5z * mg * vel[2] * sin08) + (2.0 * c5x * c5y * l5m * vel[4] * pow(cos02, 2.0)) + (
                              2.0 * cgx * cgy * mg * vel[4] * pow(cos02, 2.0)) + (c5y * l5m * l5x * vel[4] * cos02) + (
                              c4x * l4m * l4z * vel[1] * cos03) - (c4x * l4m * l4z * vel[2] * cos03) + (
                              cgy * l5x * mg * vel[4] * cos02) + (l5m * l5x * l4z * vel[1] * cos03) - (
                              l5m * l5x * l4z * vel[2] * cos03) + (l5x * l4z * mg * vel[1] * cos03) - (
                              l5x * l4z * mg * vel[2] * cos03) + (c5x * l5m * l5x * vel[4] * sin07) - (
                              c4z * l4m * l4z * vel[1] * sin08) + (c4z * l4m * l4z * vel[2] * sin08) - (
                              c5z * l5m * l4z * vel[1] * sin08) + (c5z * l5m * l4z * vel[2] * sin08) + (
                              cgx * l5x * mg * vel[4] * sin07) - (cgz * l4z * mg * vel[1] * sin08) + (
                              cgz * l4z * mg * vel[2] * sin08) - (l5m * l4z * l5z * vel[1] * sin08) + (
                              l5m * l4z * l5z * vel[2] * sin08) + (c5x * l5m * l4z * vel[1] * cos03 * cos02) - (
                              c5x * l5m * l4z * vel[2] * cos03 * cos02) + (cgx * l4z * mg * vel[1] * cos03 * cos02) - (
                              cgx * l4z * mg * vel[2] * cos03 * cos02) - (c5y * l5m * l4z * vel[1] * cos03 * sin07) + (
                              c5y * l5m * l4z * vel[2] * cos03 * sin07) - (cgy * l4z * mg * vel[1] * cos03 * sin07) + (
                              cgy * l4z * mg * vel[2] * cos03 * sin07)
        C[3, 3] = -(vel[4] * ((I5y * sin00) - (I5x * sin00) - (Igx * sin00) + (Igy * sin00) + (
                    2.0 * l5m * ((c5y * cos02) + (c5x * sin07)) * (l5x + (c5x * cos02) - (c5y * sin07))) + (
                                          2.0 * mg * ((cgy * cos02) + (cgx * sin07)) * (
                                              l5x + (cgx * cos02) - (cgy * sin07)))) / 2.0)
        C[3, 4] = (I5x * vel[1] * sin00 / 2.0) - (I5x * vel[2] * sin00 / 2.0) + (I5x * vel[3] * sin00 / 2.0) - (
                    I5y * vel[1] * sin00 / 2.0) + (I5y * vel[2] * sin00 / 2.0) - (I5y * vel[3] * sin00 / 2.0) + (
                              Igx * vel[1] * sin00 / 2.0) - (Igx * vel[2] * sin00 / 2.0) + (
                              Igx * vel[3] * sin00 / 2.0) - (Igy * vel[1] * sin00 / 2.0) + (
                              Igy * vel[2] * sin00 / 2.0) - (Igy * vel[3] * sin00 / 2.0) - (
                              I5x * vel[0] * sin01 / 4.0) - (I5x * vel[0] * sin02 / 4.0) + (
                              I5y * vel[0] * sin01 / 4.0) + (I5y * vel[0] * sin02 / 4.0) - (
                              Igx * vel[0] * sin01 / 4.0) - (Igx * vel[0] * sin02 / 4.0) + (
                              Igy * vel[0] * sin01 / 4.0) + (Igy * vel[0] * sin02 / 4.0) + (
                              I5z * vel[0] * sin06 / 2.0) + (Igz * vel[0] * sin06 / 2.0) - (
                              pow(c5x, 2.0) * l5m * vel[1] * sin00 / 2.0) + (
                              pow(c5x, 2.0) * l5m * vel[2] * sin00 / 2.0) - (
                              pow(c5x, 2.0) * l5m * vel[3] * sin00 / 2.0) + (
                              pow(c5y, 2.0) * l5m * vel[1] * sin00 / 2.0) - (
                              pow(c5y, 2.0) * l5m * vel[2] * sin00 / 2.0) + (
                              pow(c5y, 2.0) * l5m * vel[3] * sin00 / 2.0) - (
                              pow(cgx, 2.0) * mg * vel[1] * sin00 / 2.0) + (
                              pow(cgx, 2.0) * mg * vel[2] * sin00 / 2.0) - (
                              pow(cgx, 2.0) * mg * vel[3] * sin00 / 2.0) + (
                              pow(cgy, 2.0) * mg * vel[1] * sin00 / 2.0) - (
                              pow(cgy, 2.0) * mg * vel[2] * sin00 / 2.0) + (
                              pow(cgy, 2.0) * mg * vel[3] * sin00 / 2.0) + (
                              pow(c5x, 2.0) * l5m * vel[0] * sin01 / 4.0) + (
                              pow(c5x, 2.0) * l5m * vel[0] * sin02 / 4.0) - (
                              pow(c5y, 2.0) * l5m * vel[0] * sin01 / 4.0) - (
                              pow(c5y, 2.0) * l5m * vel[0] * sin02 / 4.0) + (
                              pow(cgx, 2.0) * mg * vel[0] * sin01 / 4.0) + (
                              pow(cgx, 2.0) * mg * vel[0] * sin02 / 4.0) - (
                              pow(cgy, 2.0) * mg * vel[0] * sin01 / 4.0) - (
                              pow(cgy, 2.0) * mg * vel[0] * sin02 / 4.0) + (
                              pow(c5x, 2.0) * l5m * vel[0] * sin06 / 2.0) + (
                              pow(c5y, 2.0) * l5m * vel[0] * sin06 / 2.0) + (
                              pow(cgx, 2.0) * mg * vel[0] * sin06 / 2.0) + (
                              pow(cgy, 2.0) * mg * vel[0] * sin06 / 2.0) + (c5y * l5m * l5x * vel[0] * cos01 / 2.0) + (
                              c5x * l5m * l5z * vel[0] * cos01 / 2.0) + (cgy * l5x * mg * vel[0] * cos01 / 2.0) + (
                              cgx * l5z * mg * vel[0] * cos01 / 2.0) - (c5y * c5z * l5m * vel[0] * sin03 / 2.0) - (
                              cgy * cgz * mg * vel[0] * sin03 / 2.0) + (c5x * l5m * l5x * vel[0] * sin03 / 2.0) - (
                              c5y * l5m * l5z * vel[0] * sin03 / 2.0) + (cgx * l5x * mg * vel[0] * sin03 / 2.0) - (
                              cgy * l5z * mg * vel[0] * sin03 / 2.0) - (c5x * c5y * l5m * vel[1] * cos11) + (
                              c5x * c5y * l5m * vel[2] * cos11) - (c5x * c5y * l5m * vel[3] * cos11) - (
                              cgx * cgy * mg * vel[1] * cos11) + (cgx * cgy * mg * vel[2] * cos11) - (
                              cgx * cgy * mg * vel[3] * cos11) - (c5x * c5y * l5m * vel[0] * cos08 / 2.0) + (
                              c5x * c5y * l5m * vel[0] * cos10 / 2.0) + (c5x * c5z * l5m * vel[0] * cos00 / 2.0) - (
                              cgx * cgy * mg * vel[0] * cos08 / 2.0) + (cgx * cgy * mg * vel[0] * cos10 / 2.0) + (
                              cgx * cgz * mg * vel[0] * cos00 / 2.0) - (c5y * l5m * l5x * vel[0] * cos00 / 2.0) + (
                              c5x * l5m * l5z * vel[0] * cos00 / 2.0) - (cgy * l5x * mg * vel[0] * cos00 / 2.0) + (
                              cgx * l5z * mg * vel[0] * cos00 / 2.0) + (c5y * c5z * l5m * vel[0] * sin04 / 2.0) + (
                              cgy * cgz * mg * vel[0] * sin04 / 2.0) + (c5x * l5m * l5x * vel[0] * sin04 / 2.0) + (
                              c5y * l5m * l5z * vel[0] * sin04 / 2.0) + (cgx * l5x * mg * vel[0] * sin04 / 2.0) + (
                              cgy * l5z * mg * vel[0] * sin04 / 2.0) + (c5x * c5z * l5m * vel[4] * cos02) + (
                              cgx * cgz * mg * vel[4] * cos02) - (c5y * l5m * l5x * vel[1] * cos02) + (
                              c5y * l5m * l5x * vel[2] * cos02) - (c5y * l5m * l5x * vel[3] * cos02) + (
                              c5x * l5m * l5z * vel[4] * cos02) - (cgy * l5x * mg * vel[1] * cos02) + (
                              cgy * l5x * mg * vel[2] * cos02) - (cgy * l5x * mg * vel[3] * cos02) + (
                              cgx * l5z * mg * vel[4] * cos02) - (c5y * c5z * l5m * vel[4] * sin07) - (
                              cgy * cgz * mg * vel[4] * sin07) - (c5x * l5m * l5x * vel[1] * sin07) + (
                              c5x * l5m * l5x * vel[2] * sin07) - (c5x * l5m * l5x * vel[3] * sin07) - (
                              c5y * l5m * l5z * vel[4] * sin07) - (cgx * l5x * mg * vel[1] * sin07) + (
                              cgx * l5x * mg * vel[2] * sin07) - (cgx * l5x * mg * vel[3] * sin07) - (
                              cgy * l5z * mg * vel[4] * sin07) + (c5x * c5z * l5m * vel[0] * cos01 / 2.0) + (
                              cgx * cgz * mg * vel[0] * cos01 / 2.0)
        C[4, 0] = (I5x * vel[0] * sin10 / 8.0) - (I5x * vel[0] * sin09 / 8.0) - (I5y * vel[0] * sin10 / 8.0) + (
                    I5y * vel[0] * sin09 / 8.0) + (Igx * vel[0] * sin10 / 8.0) - (Igx * vel[0] * sin09 / 8.0) - (
                              Igy * vel[0] * sin10 / 8.0) + (Igy * vel[0] * sin09 / 8.0) + (
                              I5x * vel[0] * sin00 / 4.0) - (I5y * vel[0] * sin00 / 4.0) + (
                              Igx * vel[0] * sin00 / 4.0) - (Igy * vel[0] * sin00 / 4.0) + (
                              I5x * vel[1] * sin01 / 4.0) + (I5x * vel[1] * sin02 / 4.0) - (
                              I5x * vel[2] * sin01 / 4.0) - (I5x * vel[2] * sin02 / 4.0) + (
                              I5x * vel[3] * sin01 / 4.0) + (I5x * vel[3] * sin02 / 4.0) - (
                              I5y * vel[1] * sin01 / 4.0) - (I5y * vel[1] * sin02 / 4.0) + (
                              I5y * vel[2] * sin01 / 4.0) + (I5y * vel[2] * sin02 / 4.0) - (
                              I5y * vel[3] * sin01 / 4.0) - (I5y * vel[3] * sin02 / 4.0) + (
                              Igx * vel[1] * sin01 / 4.0) + (Igx * vel[1] * sin02 / 4.0) - (
                              Igx * vel[2] * sin01 / 4.0) - (Igx * vel[2] * sin02 / 4.0) + (
                              Igx * vel[3] * sin01 / 4.0) + (Igx * vel[3] * sin02 / 4.0) - (
                              Igy * vel[1] * sin01 / 4.0) - (Igy * vel[1] * sin02 / 4.0) + (
                              Igy * vel[2] * sin01 / 4.0) + (Igy * vel[2] * sin02 / 4.0) - (
                              Igy * vel[3] * sin01 / 4.0) - (Igy * vel[3] * sin02 / 4.0) - (
                              I5z * vel[1] * sin06 / 2.0) + (I5z * vel[2] * sin06 / 2.0) - (
                              I5z * vel[3] * sin06 / 2.0) - (Igz * vel[1] * sin06 / 2.0) + (
                              Igz * vel[2] * sin06 / 2.0) - (Igz * vel[3] * sin06 / 2.0) - (
                              pow(c5x, 2.0) * l5m * vel[0] * sin10 / 8.0) + (
                              pow(c5x, 2.0) * l5m * vel[0] * sin09 / 8.0) + (
                              pow(c5y, 2.0) * l5m * vel[0] * sin10 / 8.0) - (
                              pow(c5y, 2.0) * l5m * vel[0] * sin09 / 8.0) - (
                              pow(cgx, 2.0) * mg * vel[0] * sin10 / 8.0) + (
                              pow(cgx, 2.0) * mg * vel[0] * sin09 / 8.0) + (
                              pow(cgy, 2.0) * mg * vel[0] * sin10 / 8.0) - (
                              pow(cgy, 2.0) * mg * vel[0] * sin09 / 8.0) - (
                              pow(c5x, 2.0) * l5m * vel[0] * sin00 / 4.0) + (
                              pow(c5y, 2.0) * l5m * vel[0] * sin00 / 4.0) - (
                              pow(cgx, 2.0) * mg * vel[0] * sin00 / 4.0) + (
                              pow(cgy, 2.0) * mg * vel[0] * sin00 / 4.0) - (
                              pow(c5x, 2.0) * l5m * vel[1] * sin01 / 4.0) - (
                              pow(c5x, 2.0) * l5m * vel[1] * sin02 / 4.0) + (
                              pow(c5x, 2.0) * l5m * vel[2] * sin01 / 4.0) + (
                              pow(c5x, 2.0) * l5m * vel[2] * sin02 / 4.0) - (
                              pow(c5x, 2.0) * l5m * vel[3] * sin01 / 4.0) - (
                              pow(c5x, 2.0) * l5m * vel[3] * sin02 / 4.0) + (
                              pow(c5y, 2.0) * l5m * vel[1] * sin01 / 4.0) + (
                              pow(c5y, 2.0) * l5m * vel[1] * sin02 / 4.0) - (
                              pow(c5y, 2.0) * l5m * vel[2] * sin01 / 4.0) - (
                              pow(c5y, 2.0) * l5m * vel[2] * sin02 / 4.0) + (
                              pow(c5y, 2.0) * l5m * vel[3] * sin01 / 4.0) + (
                              pow(c5y, 2.0) * l5m * vel[3] * sin02 / 4.0) - (
                              pow(cgx, 2.0) * mg * vel[1] * sin01 / 4.0) - (
                              pow(cgx, 2.0) * mg * vel[1] * sin02 / 4.0) + (
                              pow(cgx, 2.0) * mg * vel[2] * sin01 / 4.0) + (
                              pow(cgx, 2.0) * mg * vel[2] * sin02 / 4.0) - (
                              pow(cgx, 2.0) * mg * vel[3] * sin01 / 4.0) - (
                              pow(cgx, 2.0) * mg * vel[3] * sin02 / 4.0) + (
                              pow(cgy, 2.0) * mg * vel[1] * sin01 / 4.0) + (
                              pow(cgy, 2.0) * mg * vel[1] * sin02 / 4.0) - (
                              pow(cgy, 2.0) * mg * vel[2] * sin01 / 4.0) - (
                              pow(cgy, 2.0) * mg * vel[2] * sin02 / 4.0) + (
                              pow(cgy, 2.0) * mg * vel[3] * sin01 / 4.0) + (
                              pow(cgy, 2.0) * mg * vel[3] * sin02 / 4.0) - (
                              pow(c5x, 2.0) * l5m * vel[1] * sin06 / 2.0) + (
                              pow(c5x, 2.0) * l5m * vel[2] * sin06 / 2.0) - (
                              pow(c5x, 2.0) * l5m * vel[3] * sin06 / 2.0) - (
                              pow(c5y, 2.0) * l5m * vel[1] * sin06 / 2.0) + (
                              pow(c5y, 2.0) * l5m * vel[2] * sin06 / 2.0) - (
                              pow(c5y, 2.0) * l5m * vel[3] * sin06 / 2.0) - (
                              pow(cgx, 2.0) * mg * vel[1] * sin06 / 2.0) + (
                              pow(cgx, 2.0) * mg * vel[2] * sin06 / 2.0) - (
                              pow(cgx, 2.0) * mg * vel[3] * sin06 / 2.0) - (
                              pow(cgy, 2.0) * mg * vel[1] * sin06 / 2.0) + (
                              pow(cgy, 2.0) * mg * vel[2] * sin06 / 2.0) - (
                              pow(cgy, 2.0) * mg * vel[3] * sin06 / 2.0) + (c5y * l5m * l2x * vel[0] * cos01 / 2.0) - (
                              c5y * l5m * l5x * vel[1] * cos01 / 2.0) + (c5y * l5m * l5x * vel[2] * cos01 / 2.0) - (
                              c5y * l5m * l5x * vel[3] * cos01 / 2.0) - (c5x * l5m * l5z * vel[1] * cos01 / 2.0) + (
                              c5x * l5m * l5z * vel[2] * cos01 / 2.0) - (c5x * l5m * l5z * vel[3] * cos01 / 2.0) + (
                              cgy * l2x * mg * vel[0] * cos01 / 2.0) - (cgy * l5x * mg * vel[1] * cos01 / 2.0) + (
                              cgy * l5x * mg * vel[2] * cos01 / 2.0) - (cgy * l5x * mg * vel[3] * cos01 / 2.0) - (
                              cgx * l5z * mg * vel[1] * cos01 / 2.0) + (cgx * l5z * mg * vel[2] * cos01 / 2.0) - (
                              cgx * l5z * mg * vel[3] * cos01 / 2.0) + (c5y * c5z * l5m * vel[1] * sin03 / 2.0) - (
                              c5y * c5z * l5m * vel[2] * sin03 / 2.0) + (c5y * c5z * l5m * vel[3] * sin03 / 2.0) + (
                              cgy * cgz * mg * vel[1] * sin03 / 2.0) - (cgy * cgz * mg * vel[2] * sin03 / 2.0) + (
                              cgy * cgz * mg * vel[3] * sin03 / 2.0) + (c5x * l5m * l2x * vel[0] * sin03 / 2.0) - (
                              c5x * l5m * l5x * vel[1] * sin03 / 2.0) + (c5x * l5m * l5x * vel[2] * sin03 / 2.0) - (
                              c5x * l5m * l5x * vel[3] * sin03 / 2.0) + (c5y * l5m * l5z * vel[1] * sin03 / 2.0) - (
                              c5y * l5m * l5z * vel[2] * sin03 / 2.0) + (c5y * l5m * l5z * vel[3] * sin03 / 2.0) + (
                              cgx * l2x * mg * vel[0] * sin03 / 2.0) - (cgx * l5x * mg * vel[1] * sin03 / 2.0) + (
                              cgx * l5x * mg * vel[2] * sin03 / 2.0) - (cgx * l5x * mg * vel[3] * sin03 / 2.0) + (
                              cgy * l5z * mg * vel[1] * sin03 / 2.0) - (cgy * l5z * mg * vel[2] * sin03 / 2.0) + (
                              cgy * l5z * mg * vel[3] * sin03 / 2.0) + (c5x * c5y * l5m * vel[0] * cos23 / 4.0) + (
                              c5x * c5y * l5m * vel[0] * cos26 / 4.0) - (c5x * c5z * l5m * vel[0] * cos05 / 4.0) + (
                              cgx * cgy * mg * vel[0] * cos23 / 4.0) + (cgx * cgy * mg * vel[0] * cos26 / 4.0) - (
                              cgx * cgz * mg * vel[0] * cos05 / 4.0) - (c5x * l5m * l3z * vel[1] * cos29 / 2.0) + (
                              c5x * l5m * l4z * vel[0] * cos09 / 4.0) - (cgx * l3z * mg * vel[1] * cos29 / 2.0) + (
                              cgx * l4z * mg * vel[0] * cos09 / 4.0) + (c5y * l5m * l5x * vel[0] * cos05 / 4.0) - (
                              c5x * l5m * l5z * vel[0] * cos05 / 4.0) + (cgy * l5x * mg * vel[0] * cos05 / 4.0) - (
                              cgx * l5z * mg * vel[0] * cos05 / 4.0) - (c5y * c5z * l5m * vel[0] * sin12 / 4.0) - (
                              cgy * cgz * mg * vel[0] * sin12 / 4.0) - (c5y * l5m * l3z * vel[1] * sin32 / 2.0) + (
                              c5y * l5m * l4z * vel[0] * sin15 / 4.0) - (cgy * l3z * mg * vel[1] * sin32 / 2.0) + (
                              cgy * l4z * mg * vel[0] * sin15 / 4.0) - (c5x * l5m * l5x * vel[0] * sin12 / 4.0) - (
                              c5y * l5m * l5z * vel[0] * sin12 / 4.0) - (cgx * l5x * mg * vel[0] * sin12 / 4.0) - (
                              cgy * l5z * mg * vel[0] * sin12 / 4.0) - (c5x * c5y * l5m * vel[0] * cos11 / 2.0) - (
                              cgx * cgy * mg * vel[0] * cos11 / 2.0) + (c5x * c5y * l5m * vel[1] * cos08 / 2.0) - (
                              c5x * c5y * l5m * vel[1] * cos10 / 2.0) - (c5x * c5y * l5m * vel[2] * cos08 / 2.0) + (
                              c5x * c5y * l5m * vel[2] * cos10 / 2.0) + (c5x * c5y * l5m * vel[3] * cos08 / 2.0) - (
                              c5x * c5y * l5m * vel[3] * cos10 / 2.0) - (c5x * c5z * l5m * vel[1] * cos00 / 2.0) + (
                              c5x * c5z * l5m * vel[2] * cos00 / 2.0) - (c5x * c5z * l5m * vel[3] * cos00 / 2.0) + (
                              cgx * cgy * mg * vel[1] * cos08 / 2.0) - (cgx * cgy * mg * vel[1] * cos10 / 2.0) - (
                              cgx * cgy * mg * vel[2] * cos08 / 2.0) + (cgx * cgy * mg * vel[2] * cos10 / 2.0) + (
                              cgx * cgy * mg * vel[3] * cos08 / 2.0) - (cgx * cgy * mg * vel[3] * cos10 / 2.0) - (
                              cgx * cgz * mg * vel[1] * cos00 / 2.0) + (cgx * cgz * mg * vel[2] * cos00 / 2.0) - (
                              cgx * cgz * mg * vel[3] * cos00 / 2.0) + (c5y * l5m * l2x * vel[0] * cos00 / 2.0) + (
                              c5y * l5m * l5x * vel[1] * cos00 / 2.0) - (c5y * l5m * l5x * vel[2] * cos00 / 2.0) + (
                              c5y * l5m * l5x * vel[3] * cos00 / 2.0) + (c5x * l5m * l3z * vel[0] * cos25 / 4.0) + (
                              c5x * l5m * l4z * vel[0] * cos24 / 4.0) - (c5x * l5m * l5z * vel[1] * cos00 / 2.0) + (
                              c5x * l5m * l5z * vel[2] * cos00 / 2.0) - (c5x * l5m * l5z * vel[3] * cos00 / 2.0) + (
                              cgy * l2x * mg * vel[0] * cos00 / 2.0) + (cgy * l5x * mg * vel[1] * cos00 / 2.0) - (
                              cgy * l5x * mg * vel[2] * cos00 / 2.0) + (cgy * l5x * mg * vel[3] * cos00 / 2.0) + (
                              cgx * l3z * mg * vel[0] * cos25 / 4.0) + (cgx * l4z * mg * vel[0] * cos24 / 4.0) - (
                              cgx * l5z * mg * vel[1] * cos00 / 2.0) + (cgx * l5z * mg * vel[2] * cos00 / 2.0) - (
                              cgx * l5z * mg * vel[3] * cos00 / 2.0) - (c5y * c5z * l5m * vel[1] * sin04 / 2.0) + (
                              c5y * c5z * l5m * vel[2] * sin04 / 2.0) - (c5y * c5z * l5m * vel[3] * sin04 / 2.0) - (
                              cgy * cgz * mg * vel[1] * sin04 / 2.0) + (cgy * cgz * mg * vel[2] * sin04 / 2.0) - (
                              cgy * cgz * mg * vel[3] * sin04 / 2.0) - (c5x * l5m * l2x * vel[0] * sin04 / 2.0) - (
                              c5x * l5m * l5x * vel[1] * sin04 / 2.0) + (c5x * l5m * l5x * vel[2] * sin04 / 2.0) - (
                              c5x * l5m * l5x * vel[3] * sin04 / 2.0) - (c5y * l5m * l3z * vel[0] * sin26 / 4.0) - (
                              c5y * l5m * l4z * vel[0] * sin28 / 4.0) - (c5y * l5m * l5z * vel[1] * sin04 / 2.0) + (
                              c5y * l5m * l5z * vel[2] * sin04 / 2.0) - (c5y * l5m * l5z * vel[3] * sin04 / 2.0) - (
                              cgx * l2x * mg * vel[0] * sin04 / 2.0) - (cgx * l5x * mg * vel[1] * sin04 / 2.0) + (
                              cgx * l5x * mg * vel[2] * sin04 / 2.0) - (cgx * l5x * mg * vel[3] * sin04 / 2.0) - (
                              cgy * l3z * mg * vel[0] * sin26 / 4.0) - (cgy * l4z * mg * vel[0] * sin28 / 4.0) - (
                              cgy * l5z * mg * vel[1] * sin04 / 2.0) + (cgy * l5z * mg * vel[2] * sin04 / 2.0) - (
                              cgy * l5z * mg * vel[3] * sin04 / 2.0) + (c5x * l5m * l3z * vel[0] * cos14 / 4.0) - (
                              c5x * l5m * l4z * vel[1] * cos16 / 2.0) + (c5x * l5m * l4z * vel[2] * cos16 / 2.0) + (
                              cgx * l3z * mg * vel[0] * cos14 / 4.0) - (cgx * l4z * mg * vel[1] * cos16 / 2.0) + (
                              cgx * l4z * mg * vel[2] * cos16 / 2.0) - (c5y * l5m * l3z * vel[0] * sin16 / 4.0) + (
                              c5y * l5m * l4z * vel[1] * sin21 / 2.0) - (c5y * l5m * l4z * vel[2] * sin21 / 2.0) - (
                              cgy * l3z * mg * vel[0] * sin16 / 4.0) + (cgy * l4z * mg * vel[1] * sin21 / 2.0) - (
                              cgy * l4z * mg * vel[2] * sin21 / 2.0) + (c5x * c5z * l5m * vel[0] * cos06 / 4.0) + (
                              cgx * cgz * mg * vel[0] * cos06 / 4.0) - (c5x * l5m * l3z * vel[1] * cos28 / 2.0) - (
                              c5x * l5m * l4z * vel[0] * cos07 / 4.0) - (cgx * l3z * mg * vel[1] * cos28 / 2.0) - (
                              cgx * l4z * mg * vel[0] * cos07 / 4.0) + (c5y * l5m * l5x * vel[0] * cos06 / 4.0) - (
                              c5x * l5m * l3z * vel[0] * cos22 / 4.0) - (c5x * l5m * l4z * vel[0] * cos21 / 4.0) + (
                              c5x * l5m * l5z * vel[0] * cos06 / 4.0) + (cgy * l5x * mg * vel[0] * cos06 / 4.0) - (
                              cgx * l3z * mg * vel[0] * cos22 / 4.0) - (cgx * l4z * mg * vel[0] * cos21 / 4.0) + (
                              cgx * l5z * mg * vel[0] * cos06 / 4.0) - (c5y * c5z * l5m * vel[0] * sin13 / 4.0) - (
                              cgy * cgz * mg * vel[0] * sin13 / 4.0) + (c5y * l5m * l3z * vel[1] * sin31 / 2.0) + (
                              c5y * l5m * l4z * vel[0] * sin14 / 4.0) + (cgy * l3z * mg * vel[1] * sin31 / 2.0) + (
                              cgy * l4z * mg * vel[0] * sin14 / 4.0) + (c5x * l5m * l5x * vel[0] * sin13 / 4.0) - (
                              c5y * l5m * l3z * vel[0] * sin27 / 4.0) - (c5y * l5m * l4z * vel[0] * sin25 / 4.0) - (
                              c5y * l5m * l5z * vel[0] * sin13 / 4.0) + (cgx * l5x * mg * vel[0] * sin13 / 4.0) - (
                              cgy * l3z * mg * vel[0] * sin27 / 4.0) - (cgy * l4z * mg * vel[0] * sin25 / 4.0) - (
                              cgy * l5z * mg * vel[0] * sin13 / 4.0) + (c5y * l5m * l5x * vel[0] * cos02 / 2.0) + (
                              cgy * l5x * mg * vel[0] * cos02 / 2.0) - (c5x * l5m * l3z * vel[0] * cos13 / 4.0) - (
                              c5x * l5m * l4z * vel[1] * cos17 / 2.0) + (c5x * l5m * l4z * vel[2] * cos17 / 2.0) - (
                              cgx * l3z * mg * vel[0] * cos13 / 4.0) - (cgx * l4z * mg * vel[1] * cos17 / 2.0) + (
                              cgx * l4z * mg * vel[2] * cos17 / 2.0) + (c5x * l5m * l5x * vel[0] * sin07 / 2.0) + (
                              cgx * l5x * mg * vel[0] * sin07 / 2.0) - (c5x * c5z * l5m * vel[1] * cos01 / 2.0) + (
                              c5x * c5z * l5m * vel[2] * cos01 / 2.0) - (c5x * c5z * l5m * vel[3] * cos01 / 2.0) - (
                              cgx * cgz * mg * vel[1] * cos01 / 2.0) + (cgx * cgz * mg * vel[2] * cos01 / 2.0) - (
                              cgx * cgz * mg * vel[3] * cos01 / 2.0) - (c5y * l5m * l3z * vel[0] * sin17 / 4.0) - (
                              c5y * l5m * l4z * vel[1] * sin22 / 2.0) + (c5y * l5m * l4z * vel[2] * sin22 / 2.0) - (
                              cgy * l3z * mg * vel[0] * sin17 / 4.0) - (cgy * l4z * mg * vel[1] * sin22 / 2.0) + (
                              cgy * l4z * mg * vel[2] * sin22 / 2.0)
        C[4, 1] = (I5x * vel[2] * sin00 / 2.0) - (I5x * vel[1] * sin00 / 2.0) - (I5x * vel[3] * sin00 / 2.0) + (
                    I5y * vel[1] * sin00 / 2.0) - (I5y * vel[2] * sin00 / 2.0) + (I5y * vel[3] * sin00 / 2.0) - (
                              Igx * vel[1] * sin00 / 2.0) + (Igx * vel[2] * sin00 / 2.0) - (
                              Igx * vel[3] * sin00 / 2.0) + (Igy * vel[1] * sin00 / 2.0) - (
                              Igy * vel[2] * sin00 / 2.0) + (Igy * vel[3] * sin00 / 2.0) + (
                              I5x * vel[0] * sin01 / 4.0) + (I5x * vel[0] * sin02 / 4.0) - (
                              I5y * vel[0] * sin01 / 4.0) - (I5y * vel[0] * sin02 / 4.0) + (
                              Igx * vel[0] * sin01 / 4.0) + (Igx * vel[0] * sin02 / 4.0) - (
                              Igy * vel[0] * sin01 / 4.0) - (Igy * vel[0] * sin02 / 4.0) - (
                              I5z * vel[0] * sin06 / 2.0) - (Igz * vel[0] * sin06 / 2.0) + (
                              pow(c5x, 2.0) * l5m * vel[1] * sin00 / 2.0) - (
                              pow(c5x, 2.0) * l5m * vel[2] * sin00 / 2.0) + (
                              pow(c5x, 2.0) * l5m * vel[3] * sin00 / 2.0) - (
                              pow(c5y, 2.0) * l5m * vel[1] * sin00 / 2.0) + (
                              pow(c5y, 2.0) * l5m * vel[2] * sin00 / 2.0) - (
                              pow(c5y, 2.0) * l5m * vel[3] * sin00 / 2.0) + (
                              pow(cgx, 2.0) * mg * vel[1] * sin00 / 2.0) - (
                              pow(cgx, 2.0) * mg * vel[2] * sin00 / 2.0) + (
                              pow(cgx, 2.0) * mg * vel[3] * sin00 / 2.0) - (
                              pow(cgy, 2.0) * mg * vel[1] * sin00 / 2.0) + (
                              pow(cgy, 2.0) * mg * vel[2] * sin00 / 2.0) - (
                              pow(cgy, 2.0) * mg * vel[3] * sin00 / 2.0) - (
                              pow(c5x, 2.0) * l5m * vel[0] * sin01 / 4.0) - (
                              pow(c5x, 2.0) * l5m * vel[0] * sin02 / 4.0) + (
                              pow(c5y, 2.0) * l5m * vel[0] * sin01 / 4.0) + (
                              pow(c5y, 2.0) * l5m * vel[0] * sin02 / 4.0) - (
                              pow(cgx, 2.0) * mg * vel[0] * sin01 / 4.0) - (
                              pow(cgx, 2.0) * mg * vel[0] * sin02 / 4.0) + (
                              pow(cgy, 2.0) * mg * vel[0] * sin01 / 4.0) + (
                              pow(cgy, 2.0) * mg * vel[0] * sin02 / 4.0) - (
                              pow(c5x, 2.0) * l5m * vel[0] * sin06 / 2.0) - (
                              pow(c5y, 2.0) * l5m * vel[0] * sin06 / 2.0) - (
                              pow(cgx, 2.0) * mg * vel[0] * sin06 / 2.0) - (
                              pow(cgy, 2.0) * mg * vel[0] * sin06 / 2.0) - (c5y * l5m * l5x * vel[0] * cos01 / 2.0) - (
                              c5x * l5m * l5z * vel[0] * cos01 / 2.0) - (cgy * l5x * mg * vel[0] * cos01 / 2.0) - (
                              cgx * l5z * mg * vel[0] * cos01 / 2.0) + (c5y * c5z * l5m * vel[0] * sin03 / 2.0) + (
                              cgy * cgz * mg * vel[0] * sin03 / 2.0) - (c5x * l5m * l5x * vel[0] * sin03 / 2.0) + (
                              c5y * l5m * l5z * vel[0] * sin03 / 2.0) - (cgx * l5x * mg * vel[0] * sin03 / 2.0) + (
                              cgy * l5z * mg * vel[0] * sin03 / 2.0) - (c5x * l5m * l3z * vel[0] * cos29 / 2.0) + (
                              c5x * l5m * l4z * vel[1] * cos09 / 2.0) - (c5x * l5m * l4z * vel[2] * cos09 / 2.0) - (
                              cgx * l3z * mg * vel[0] * cos29 / 2.0) + (cgx * l4z * mg * vel[1] * cos09 / 2.0) - (
                              cgx * l4z * mg * vel[2] * cos09 / 2.0) - (c5y * l5m * l3z * vel[0] * sin32 / 2.0) + (
                              c5y * l5m * l4z * vel[1] * sin15 / 2.0) - (c5y * l5m * l4z * vel[2] * sin15 / 2.0) - (
                              cgy * l3z * mg * vel[0] * sin32 / 2.0) + (cgy * l4z * mg * vel[1] * sin15 / 2.0) - (
                              cgy * l4z * mg * vel[2] * sin15 / 2.0) + (c5x * c5y * l5m * vel[1] * cos11) - (
                              c5x * c5y * l5m * vel[2] * cos11) + (c5x * c5y * l5m * vel[3] * cos11) + (
                              cgx * cgy * mg * vel[1] * cos11) - (cgx * cgy * mg * vel[2] * cos11) + (
                              cgx * cgy * mg * vel[3] * cos11) + (c5x * c5y * l5m * vel[0] * cos08 / 2.0) - (
                              c5x * c5y * l5m * vel[0] * cos10 / 2.0) - (c5x * c5z * l5m * vel[0] * cos00 / 2.0) + (
                              cgx * cgy * mg * vel[0] * cos08 / 2.0) - (cgx * cgy * mg * vel[0] * cos10 / 2.0) - (
                              cgx * cgz * mg * vel[0] * cos00 / 2.0) + (c5y * l5m * l5x * vel[0] * cos00 / 2.0) - (
                              c5x * l5m * l5z * vel[0] * cos00 / 2.0) + (cgy * l5x * mg * vel[0] * cos00 / 2.0) - (
                              cgx * l5z * mg * vel[0] * cos00 / 2.0) - (c5y * c5z * l5m * vel[0] * sin04 / 2.0) - (
                              cgy * cgz * mg * vel[0] * sin04 / 2.0) - (c5x * l5m * l5x * vel[0] * sin04 / 2.0) - (
                              c5y * l5m * l5z * vel[0] * sin04 / 2.0) - (cgx * l5x * mg * vel[0] * sin04 / 2.0) - (
                              cgy * l5z * mg * vel[0] * sin04 / 2.0) - (c5x * l5m * l4z * vel[0] * cos16 / 2.0) + (
                              c5x * l5m * l3z * vel[1] * cos14 / 2.0) - (cgx * l4z * mg * vel[0] * cos16 / 2.0) + (
                              cgx * l3z * mg * vel[1] * cos14 / 2.0) + (c5y * l5m * l4z * vel[0] * sin21 / 2.0) - (
                              c5y * l5m * l3z * vel[1] * sin16 / 2.0) + (cgy * l4z * mg * vel[0] * sin21 / 2.0) - (
                              cgy * l3z * mg * vel[1] * sin16 / 2.0) - (c5x * l5m * l3z * vel[0] * cos28 / 2.0) - (
                              c5x * l5m * l4z * vel[1] * cos07 / 2.0) + (c5x * l5m * l4z * vel[2] * cos07 / 2.0) - (
                              cgx * l3z * mg * vel[0] * cos28 / 2.0) - (cgx * l4z * mg * vel[1] * cos07 / 2.0) + (
                              cgx * l4z * mg * vel[2] * cos07 / 2.0) + (c5y * l5m * l3z * vel[0] * sin31 / 2.0) + (
                              c5y * l5m * l4z * vel[1] * sin14 / 2.0) - (c5y * l5m * l4z * vel[2] * sin14 / 2.0) + (
                              cgy * l3z * mg * vel[0] * sin31 / 2.0) + (cgy * l4z * mg * vel[1] * sin14 / 2.0) - (
                              cgy * l4z * mg * vel[2] * sin14 / 2.0) + (c5y * l5m * l5x * vel[1] * cos02) - (
                              c5y * l5m * l5x * vel[2] * cos02) + (c5y * l5m * l5x * vel[3] * cos02) + (
                              cgy * l5x * mg * vel[1] * cos02) - (cgy * l5x * mg * vel[2] * cos02) + (
                              cgy * l5x * mg * vel[3] * cos02) - (c5x * l5m * l4z * vel[0] * cos17 / 2.0) - (
                              c5x * l5m * l3z * vel[1] * cos13 / 2.0) - (cgx * l4z * mg * vel[0] * cos17 / 2.0) - (
                              cgx * l3z * mg * vel[1] * cos13 / 2.0) + (c5x * l5m * l5x * vel[1] * sin07) - (
                              c5x * l5m * l5x * vel[2] * sin07) + (c5x * l5m * l5x * vel[3] * sin07) + (
                              cgx * l5x * mg * vel[1] * sin07) - (cgx * l5x * mg * vel[2] * sin07) + (
                              cgx * l5x * mg * vel[3] * sin07) - (c5x * c5z * l5m * vel[0] * cos01 / 2.0) - (
                              cgx * cgz * mg * vel[0] * cos01 / 2.0) - (c5y * l5m * l4z * vel[0] * sin22 / 2.0) - (
                              c5y * l5m * l3z * vel[1] * sin17 / 2.0) - (cgy * l4z * mg * vel[0] * sin22 / 2.0) - (
                              cgy * l3z * mg * vel[1] * sin17 / 2.0)
        C[4, 2] = (I5x * vel[1] * sin00 / 2.0) - (I5x * vel[2] * sin00 / 2.0) + (I5x * vel[3] * sin00 / 2.0) - (
                    I5y * vel[1] * sin00 / 2.0) + (I5y * vel[2] * sin00 / 2.0) - (I5y * vel[3] * sin00 / 2.0) + (
                              Igx * vel[1] * sin00 / 2.0) - (Igx * vel[2] * sin00 / 2.0) + (
                              Igx * vel[3] * sin00 / 2.0) - (Igy * vel[1] * sin00 / 2.0) + (
                              Igy * vel[2] * sin00 / 2.0) - (Igy * vel[3] * sin00 / 2.0) - (
                              I5x * vel[0] * sin01 / 4.0) - (I5x * vel[0] * sin02 / 4.0) + (
                              I5y * vel[0] * sin01 / 4.0) + (I5y * vel[0] * sin02 / 4.0) - (
                              Igx * vel[0] * sin01 / 4.0) - (Igx * vel[0] * sin02 / 4.0) + (
                              Igy * vel[0] * sin01 / 4.0) + (Igy * vel[0] * sin02 / 4.0) + (
                              I5z * vel[0] * sin06 / 2.0) + (Igz * vel[0] * sin06 / 2.0) - (
                              pow(c5x, 2.0) * l5m * vel[1] * sin00 / 2.0) + (
                              pow(c5x, 2.0) * l5m * vel[2] * sin00 / 2.0) - (
                              pow(c5x, 2.0) * l5m * vel[3] * sin00 / 2.0) + (
                              pow(c5y, 2.0) * l5m * vel[1] * sin00 / 2.0) - (
                              pow(c5y, 2.0) * l5m * vel[2] * sin00 / 2.0) + (
                              pow(c5y, 2.0) * l5m * vel[3] * sin00 / 2.0) - (
                              pow(cgx, 2.0) * mg * vel[1] * sin00 / 2.0) + (
                              pow(cgx, 2.0) * mg * vel[2] * sin00 / 2.0) - (
                              pow(cgx, 2.0) * mg * vel[3] * sin00 / 2.0) + (
                              pow(cgy, 2.0) * mg * vel[1] * sin00 / 2.0) - (
                              pow(cgy, 2.0) * mg * vel[2] * sin00 / 2.0) + (
                              pow(cgy, 2.0) * mg * vel[3] * sin00 / 2.0) + (
                              pow(c5x, 2.0) * l5m * vel[0] * sin01 / 4.0) + (
                              pow(c5x, 2.0) * l5m * vel[0] * sin02 / 4.0) - (
                              pow(c5y, 2.0) * l5m * vel[0] * sin01 / 4.0) - (
                              pow(c5y, 2.0) * l5m * vel[0] * sin02 / 4.0) + (
                              pow(cgx, 2.0) * mg * vel[0] * sin01 / 4.0) + (
                              pow(cgx, 2.0) * mg * vel[0] * sin02 / 4.0) - (
                              pow(cgy, 2.0) * mg * vel[0] * sin01 / 4.0) - (
                              pow(cgy, 2.0) * mg * vel[0] * sin02 / 4.0) + (
                              pow(c5x, 2.0) * l5m * vel[0] * sin06 / 2.0) + (
                              pow(c5y, 2.0) * l5m * vel[0] * sin06 / 2.0) + (
                              pow(cgx, 2.0) * mg * vel[0] * sin06 / 2.0) + (
                              pow(cgy, 2.0) * mg * vel[0] * sin06 / 2.0) + (c5y * l5m * l5x * vel[0] * cos01 / 2.0) + (
                              c5x * l5m * l5z * vel[0] * cos01 / 2.0) + (cgy * l5x * mg * vel[0] * cos01 / 2.0) + (
                              cgx * l5z * mg * vel[0] * cos01 / 2.0) - (c5y * c5z * l5m * vel[0] * sin03 / 2.0) - (
                              cgy * cgz * mg * vel[0] * sin03 / 2.0) + (c5x * l5m * l5x * vel[0] * sin03 / 2.0) - (
                              c5y * l5m * l5z * vel[0] * sin03 / 2.0) + (cgx * l5x * mg * vel[0] * sin03 / 2.0) - (
                              cgy * l5z * mg * vel[0] * sin03 / 2.0) - (c5x * l5m * l4z * vel[1] * cos09 / 2.0) + (
                              c5x * l5m * l4z * vel[2] * cos09 / 2.0) - (cgx * l4z * mg * vel[1] * cos09 / 2.0) + (
                              cgx * l4z * mg * vel[2] * cos09 / 2.0) - (c5y * l5m * l4z * vel[1] * sin15 / 2.0) + (
                              c5y * l5m * l4z * vel[2] * sin15 / 2.0) - (cgy * l4z * mg * vel[1] * sin15 / 2.0) + (
                              cgy * l4z * mg * vel[2] * sin15 / 2.0) - (c5x * c5y * l5m * vel[1] * cos11) + (
                              c5x * c5y * l5m * vel[2] * cos11) - (c5x * c5y * l5m * vel[3] * cos11) - (
                              cgx * cgy * mg * vel[1] * cos11) + (cgx * cgy * mg * vel[2] * cos11) - (
                              cgx * cgy * mg * vel[3] * cos11) - (c5x * c5y * l5m * vel[0] * cos08 / 2.0) + (
                              c5x * c5y * l5m * vel[0] * cos10 / 2.0) + (c5x * c5z * l5m * vel[0] * cos00 / 2.0) - (
                              cgx * cgy * mg * vel[0] * cos08 / 2.0) + (cgx * cgy * mg * vel[0] * cos10 / 2.0) + (
                              cgx * cgz * mg * vel[0] * cos00 / 2.0) - (c5y * l5m * l5x * vel[0] * cos00 / 2.0) + (
                              c5x * l5m * l5z * vel[0] * cos00 / 2.0) - (cgy * l5x * mg * vel[0] * cos00 / 2.0) + (
                              cgx * l5z * mg * vel[0] * cos00 / 2.0) + (c5y * c5z * l5m * vel[0] * sin04 / 2.0) + (
                              cgy * cgz * mg * vel[0] * sin04 / 2.0) + (c5x * l5m * l5x * vel[0] * sin04 / 2.0) + (
                              c5y * l5m * l5z * vel[0] * sin04 / 2.0) + (cgx * l5x * mg * vel[0] * sin04 / 2.0) + (
                              cgy * l5z * mg * vel[0] * sin04 / 2.0) + (c5x * l5m * l4z * vel[0] * cos16 / 2.0) + (
                              cgx * l4z * mg * vel[0] * cos16 / 2.0) - (c5y * l5m * l4z * vel[0] * sin21 / 2.0) - (
                              cgy * l4z * mg * vel[0] * sin21 / 2.0) + (c5x * l5m * l4z * vel[1] * cos07 / 2.0) - (
                              c5x * l5m * l4z * vel[2] * cos07 / 2.0) + (cgx * l4z * mg * vel[1] * cos07 / 2.0) - (
                              cgx * l4z * mg * vel[2] * cos07 / 2.0) - (c5y * l5m * l4z * vel[1] * sin14 / 2.0) + (
                              c5y * l5m * l4z * vel[2] * sin14 / 2.0) - (cgy * l4z * mg * vel[1] * sin14 / 2.0) + (
                              cgy * l4z * mg * vel[2] * sin14 / 2.0) - (c5y * l5m * l5x * vel[1] * cos02) + (
                              c5y * l5m * l5x * vel[2] * cos02) - (c5y * l5m * l5x * vel[3] * cos02) - (
                              cgy * l5x * mg * vel[1] * cos02) + (cgy * l5x * mg * vel[2] * cos02) - (
                              cgy * l5x * mg * vel[3] * cos02) + (c5x * l5m * l4z * vel[0] * cos17 / 2.0) + (
                              cgx * l4z * mg * vel[0] * cos17 / 2.0) - (c5x * l5m * l5x * vel[1] * sin07) + (
                              c5x * l5m * l5x * vel[2] * sin07) - (c5x * l5m * l5x * vel[3] * sin07) - (
                              cgx * l5x * mg * vel[1] * sin07) + (cgx * l5x * mg * vel[2] * sin07) - (
                              cgx * l5x * mg * vel[3] * sin07) + (c5x * c5z * l5m * vel[0] * cos01 / 2.0) + (
                              cgx * cgz * mg * vel[0] * cos01 / 2.0) + (c5y * l5m * l4z * vel[0] * sin22 / 2.0) + (
                              cgy * l4z * mg * vel[0] * sin22 / 2.0)
        C[4, 3] = (I5x * vel[2] * sin00 / 2.0) - (I5x * vel[1] * sin00 / 2.0) - (I5x * vel[3] * sin00 / 2.0) + (
                    I5y * vel[1] * sin00 / 2.0) - (I5y * vel[2] * sin00 / 2.0) + (I5y * vel[3] * sin00 / 2.0) - (
                              Igx * vel[1] * sin00 / 2.0) + (Igx * vel[2] * sin00 / 2.0) - (
                              Igx * vel[3] * sin00 / 2.0) + (Igy * vel[1] * sin00 / 2.0) - (
                              Igy * vel[2] * sin00 / 2.0) + (Igy * vel[3] * sin00 / 2.0) + (
                              I5x * vel[0] * sin01 / 4.0) + (I5x * vel[0] * sin02 / 4.0) - (
                              I5y * vel[0] * sin01 / 4.0) - (I5y * vel[0] * sin02 / 4.0) + (
                              Igx * vel[0] * sin01 / 4.0) + (Igx * vel[0] * sin02 / 4.0) - (
                              Igy * vel[0] * sin01 / 4.0) - (Igy * vel[0] * sin02 / 4.0) - (
                              I5z * vel[0] * sin06 / 2.0) - (Igz * vel[0] * sin06 / 2.0) + (
                              pow(c5x, 2.0) * l5m * vel[1] * sin00 / 2.0) - (
                              pow(c5x, 2.0) * l5m * vel[2] * sin00 / 2.0) + (
                              pow(c5x, 2.0) * l5m * vel[3] * sin00 / 2.0) - (
                              pow(c5y, 2.0) * l5m * vel[1] * sin00 / 2.0) + (
                              pow(c5y, 2.0) * l5m * vel[2] * sin00 / 2.0) - (
                              pow(c5y, 2.0) * l5m * vel[3] * sin00 / 2.0) + (
                              pow(cgx, 2.0) * mg * vel[1] * sin00 / 2.0) - (
                              pow(cgx, 2.0) * mg * vel[2] * sin00 / 2.0) + (
                              pow(cgx, 2.0) * mg * vel[3] * sin00 / 2.0) - (
                              pow(cgy, 2.0) * mg * vel[1] * sin00 / 2.0) + (
                              pow(cgy, 2.0) * mg * vel[2] * sin00 / 2.0) - (
                              pow(cgy, 2.0) * mg * vel[3] * sin00 / 2.0) - (
                              pow(c5x, 2.0) * l5m * vel[0] * sin01 / 4.0) - (
                              pow(c5x, 2.0) * l5m * vel[0] * sin02 / 4.0) + (
                              pow(c5y, 2.0) * l5m * vel[0] * sin01 / 4.0) + (
                              pow(c5y, 2.0) * l5m * vel[0] * sin02 / 4.0) - (
                              pow(cgx, 2.0) * mg * vel[0] * sin01 / 4.0) - (
                              pow(cgx, 2.0) * mg * vel[0] * sin02 / 4.0) + (
                              pow(cgy, 2.0) * mg * vel[0] * sin01 / 4.0) + (
                              pow(cgy, 2.0) * mg * vel[0] * sin02 / 4.0) - (
                              pow(c5x, 2.0) * l5m * vel[0] * sin06 / 2.0) - (
                              pow(c5y, 2.0) * l5m * vel[0] * sin06 / 2.0) - (
                              pow(cgx, 2.0) * mg * vel[0] * sin06 / 2.0) - (
                              pow(cgy, 2.0) * mg * vel[0] * sin06 / 2.0) - (c5y * l5m * l5x * vel[0] * cos01 / 2.0) - (
                              c5x * l5m * l5z * vel[0] * cos01 / 2.0) - (cgy * l5x * mg * vel[0] * cos01 / 2.0) - (
                              cgx * l5z * mg * vel[0] * cos01 / 2.0) + (c5y * c5z * l5m * vel[0] * sin03 / 2.0) + (
                              cgy * cgz * mg * vel[0] * sin03 / 2.0) - (c5x * l5m * l5x * vel[0] * sin03 / 2.0) + (
                              c5y * l5m * l5z * vel[0] * sin03 / 2.0) - (cgx * l5x * mg * vel[0] * sin03 / 2.0) + (
                              cgy * l5z * mg * vel[0] * sin03 / 2.0) + (c5x * c5y * l5m * vel[1] * cos11) - (
                              c5x * c5y * l5m * vel[2] * cos11) + (c5x * c5y * l5m * vel[3] * cos11) + (
                              cgx * cgy * mg * vel[1] * cos11) - (cgx * cgy * mg * vel[2] * cos11) + (
                              cgx * cgy * mg * vel[3] * cos11) + (c5x * c5y * l5m * vel[0] * cos08 / 2.0) - (
                              c5x * c5y * l5m * vel[0] * cos10 / 2.0) - (c5x * c5z * l5m * vel[0] * cos00 / 2.0) + (
                              cgx * cgy * mg * vel[0] * cos08 / 2.0) - (cgx * cgy * mg * vel[0] * cos10 / 2.0) - (
                              cgx * cgz * mg * vel[0] * cos00 / 2.0) + (c5y * l5m * l5x * vel[0] * cos00 / 2.0) - (
                              c5x * l5m * l5z * vel[0] * cos00 / 2.0) + (cgy * l5x * mg * vel[0] * cos00 / 2.0) - (
                              cgx * l5z * mg * vel[0] * cos00 / 2.0) - (c5y * c5z * l5m * vel[0] * sin04 / 2.0) - (
                              cgy * cgz * mg * vel[0] * sin04 / 2.0) - (c5x * l5m * l5x * vel[0] * sin04 / 2.0) - (
                              c5y * l5m * l5z * vel[0] * sin04 / 2.0) - (cgx * l5x * mg * vel[0] * sin04 / 2.0) - (
                              cgy * l5z * mg * vel[0] * sin04 / 2.0) + (c5y * l5m * l5x * vel[1] * cos02) - (
                              c5y * l5m * l5x * vel[2] * cos02) + (c5y * l5m * l5x * vel[3] * cos02) + (
                              cgy * l5x * mg * vel[1] * cos02) - (cgy * l5x * mg * vel[2] * cos02) + (
                              cgy * l5x * mg * vel[3] * cos02) + (c5x * l5m * l5x * vel[1] * sin07) - (
                              c5x * l5m * l5x * vel[2] * sin07) + (c5x * l5m * l5x * vel[3] * sin07) + (
                              cgx * l5x * mg * vel[1] * sin07) - (cgx * l5x * mg * vel[2] * sin07) + (
                              cgx * l5x * mg * vel[3] * sin07) - (c5x * c5z * l5m * vel[0] * cos01 / 2.0) - (
                              cgx * cgz * mg * vel[0] * cos01 / 2.0)
        C[4, 4] = 0.0


    def N_vector(self):
        pos = self.pos
        N = self.N

        lbx = 0;
        lbz = 0;
        l1x = 0.024;
        l1z = 0.115;
        l2x = 0.033;
        l2z = 0;
        l3x = 0;
        l3z = 0.155;
        l4x = 0;
        l4z = 0.135;
        l5x = 0;
        l5z = 0.1136;
        lfx = 0;
        lfz = 0.05716;

        l1m = 2.351;
        l2m = 1.318;
        l3m = 0.821;
        l4m = 0.769;
        l5m = 0.687;
        mg = 0.219;

        c1x = 0.01516;
        c1y = 0.00359;
        c1z = -0.03105;
        c2x = -0.01903;
        c2y = 0.0150;
        c2z = 0.11397;
        c3x = 0.00013;
        c3y = 0.02022;
        c3z = 0.10441;
        c4x = 0.00015;
        c4y = -0.02464;
        c4z = 0.05353;
        c5x = 0;
        c5y = 0.0012;
        c5z = 0.01648;
        cgx = 0;
        cgy = 0;
        cgz = 0.0289;

        I1x = 0.0029525;
        I1y = 0.0060091;
        I1z = 0.0058821;
        I2x = 0.0031145;
        I2y = 0.0005843;
        I2z = 0.0031631;
        I3x = 0.00172767;
        I3y = 0.00041967;
        I3z = 0.0018468;
        I4x = 0.0006764;
        I4y = 0.0010573;
        I4z = 0.0006610;
        I5x = 0.0001934;
        I5y = 0.0001602;
        I5z = 0.0000689;
        Igx = 0.0002324;
        Igy = 0.0003629;
        Igz = 0.0002067;

        sin00 = sin(pos[1] - pos[2] + pos[3]) 
        cos00 = cos(pos[1] - pos[2] + pos[3])
        sin01 = sin(pos[1] - pos[2])
        cos01 = cos(pos[1] - pos[2] + pos[3] - pos[4])
        cos02 = cos(pos[1] - pos[2] + pos[3] + pos[4])
        sin02 = sin(pos[1] - pos[2] + pos[3] - pos[4])
        sin03 = sin(pos[1] - pos[2] + pos[3] + pos[4])
        sin04 = sin(pos[1])
        N[0] = 0.0
        N[1] = (981.0 * c2x * l2m * cos(pos[1]) / 100.0) - (981.0 * c2z * l2m * sin04 / 100.0) - (
                    981.0 * l3m * l3z * sin04 / 100.0) - (981.0 * l4m * l3z * sin04 / 100.0) - (
                              981.0 * l5m * l3z * sin04 / 100.0) - (981.0 * l3z * mg * sin04 / 100.0) + (
                              981.0 * c5x * l5m * cos02 / 200.0) + (981.0 * cgx * mg * cos02 / 200.0) - (
                              981.0 * c5y * l5m * sin03 / 200.0) - (981.0 * cgy * mg * sin03 / 200.0) + (
                              981.0 * c3x * l3m * cos(pos[1] - pos[2]) / 100.0) - (
                              981.0 * c3z * l3m * sin01 / 100.0) - (981.0 * l4m * l4z * sin01 / 100.0) - (
                              981.0 * l5m * l4z * sin01 / 100.0) - (981.0 * l4z * mg * sin01 / 100.0) + (
                              981.0 * c5x * l5m * cos01 / 200.0) + (981.0 * cgx * mg * cos01 / 200.0) + (
                              981.0 * c5y * l5m * sin02 / 200.0) + (981.0 * cgy * mg * sin02 / 200.0) + (
                              981.0 * c4x * l4m * cos00 / 100.0) + (981.0 * l5m * l5x * cos00 / 100.0) + (
                              981.0 * l5x * mg * cos00 / 100.0) - (981.0 * c4z * l4m * sin00 / 100.0) - (
                              981.0 * c5z * l5m * sin00 / 100.0) - (981.0 * cgz * mg * sin00 / 100.0) - (
                              981.0 * l5m * l5z * sin00 / 100.0) - (981.0 * l5z * mg * sin00 / 100.0);
        N[2] = (981.0 * c5y * l5m * sin03 / 200.0) - (981.0 * cgx * mg * cos02 / 200.0) - (
                    981.0 * c5x * l5m * cos02 / 200.0) + (981.0 * cgy * mg * sin03 / 200.0) - (
                              981.0 * c3x * l3m * cos(pos[1] - pos[2]) / 100.0) + (
                              981.0 * c3z * l3m * sin01 / 100.0) + (981.0 * l4m * l4z * sin01 / 100.0) + (
                              981.0 * l5m * l4z * sin01 / 100.0) + (981.0 * l4z * mg * sin01 / 100.0) - (
                              981.0 * c5x * l5m * cos01 / 200.0) - (981.0 * cgx * mg * cos01 / 200.0) - (
                              981.0 * c5y * l5m * sin02 / 200.0) - (981.0 * cgy * mg * sin02 / 200.0) - (
                              981.0 * c4x * l4m * cos00 / 100.0) - (981.0 * l5m * l5x * cos00 / 100.0) - (
                              981.0 * l5x * mg * cos00 / 100.0) + (981.0 * c4z * l4m * sin00 / 100.0) + (
                              981.0 * c5z * l5m * sin00 / 100.0) + (981.0 * cgz * mg * sin00 / 100.0) + (
                              981.0 * l5m * l5z * sin00 / 100.0) + (981.0 * l5z * mg * sin00 / 100.0);
        N[3] = (981.0 * c5x * l5m * cos02 / 200.0) + (981.0 * cgx * mg * cos02 / 200.0) - (
                    981.0 * c5y * l5m * sin03 / 200.0) - (981.0 * cgy * mg * sin03 / 200.0) + (
                              981.0 * c5x * l5m * cos01 / 200.0) + (981.0 * cgx * mg * cos01 / 200.0) + (
                              981.0 * c5y * l5m * sin02 / 200.0) + (981.0 * cgy * mg * sin02 / 200.0) + (
                              981.0 * c4x * l4m * cos00 / 100.0) + (981.0 * l5m * l5x * cos00 / 100.0) + (
                              981.0 * l5x * mg * cos00 / 100.0) - (981.0 * c4z * l4m * sin00 / 100.0) - (
                              981.0 * c5z * l5m * sin00 / 100.0) - (981.0 * cgz * mg * sin00 / 100.0) - (
                              981.0 * l5m * l5z * sin00 / 100.0) - (981.0 * l5z * mg * sin00 / 100.0);
        N[4] = -(981.0 * sin00 * (
                    (c5y * l5m * cos(pos[4])) + (cgy * mg * cos(pos[4])) + (c5x * l5m * sin(pos[4])) + (
                        cgx * mg * sin(pos[4]))) / 100.0)


if __name__ == '__main__':

    dy = dynamic(5)
    dy.pos = np.ones(5)
    dy.vel = np.ones(5)
    dy.C_matrix()
    dy.N_vector()
    print("C_matrix is {}".format(dy.N))

