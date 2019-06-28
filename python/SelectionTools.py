from __future__ import print_function
from __future__ import division

from ROOT import TLorentzVector



##
##  HELPERS
##


# SIGN

def sign(x):
    if x > 0:
        return 1
    elif x < 0:
        return -1
    elif x == 0:
        return 0
    else:
        return x


# Returns (vector) sum of TLorentzVectors in list
def GetSum(p4_list):
    p4 = TLorentzVector()
    for p4_ in p4_list:
        p4 += p4_
    return p4




##
##  GET_MASS_PAIRS
##

# Returns high- and low-mass pairs for same-flavor final state
# Takes TLorentzVectors (list) and charges (dict)
# Pairing is determined by maximizing mass difference


def GetMassPairs(p4_list, q_dict):

    ##  INITIALIZE  ##

    # Assign highest-p lepton as lep1, which is always in pair 1
    p4_list.sort(key=lambda fourvec: fourvec.P(), reverse=True)

    lep1p4, lep1q = p4_list[0], q_dict[p4_list[0]]
    p1_list, p2_list = [lep1p4], []

    config = 1



    ##  CONFIG 1  ##

    # OS leptons with largest momenta are Z1

    # Find subleading lepton of pair 1
    for lep in p4_list:
        if q_dict[lep] != lep1q:
            p1_list.append(lep)
            break


    # Assemble pair 2
    for lep in p4_list:
        if lep not in p1_list:
            p2_list.append(lep)


    # Find difference between masses of pairs and assign them
    pair1, pair2 = GetSum(p1_list), GetSum(p2_list)
    mass_diff_1 = pair1.M() - pair2.M()

    if mass_diff_1 > 0:
        p_pair, k_pair = pair1, pair2
        p_list, k_list = p1_list[:], p2_list[:]
    else:
        k_pair, p_pair = pair1, pair2
        k_list, p_list = p1_list[:], p2_list[:]



    ##  CONFIG 2  ##

    # Swap leptons that are not SS as lep1
    # (i.e. look for pair 2 lepton with SS as subleading lep of pair 1)
    for i in range(len(p2_list)):
        if q_dict[p2_list[i]] == q_dict[p1_list[1]]:
            p2_list[i], p1_list[1] = p1_list[1], p2_list[i]
            break


    # Reassign pair momenta
    pair1, pair2 = GetSum(p1_list), GetSum(p2_list)
    mass_diff_2 = pair1.M() - pair2.M()


    # If this config has larger mass difference, reassign returned values
    if abs(mass_diff_2) > abs(mass_diff_1):
        config = 2

        if mass_diff_2 > 0:
            p_pair, k_pair = pair1, pair2
            p_list, k_list = p1_list[:], p2_list[:]
        else:
            k_pair, p_pair = pair1, pair2
            k_list, p_list = p1_list[:], p2_list[:]



    ##  RETURN ##
    return p_pair, k_pair, p_list, k_list
