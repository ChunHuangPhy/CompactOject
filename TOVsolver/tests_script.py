import main
import EoS_import
import constant
import speed_of_sound
import TOV_solver

def test_constant():
    assert constant.c == 3e10
    assert constant.Msun == 1.989e33
    assert constant.dyncm2_to_MeVfm3 == 1./(1.6022e33)
    assert constant.gcm3_to_MeVfm3 == 1./(1.7827e12)
    assert constant.oneoverfm_MeV == 197.33
    assert constant.c == type(float)
    assert constant.Msun == type(float)
    assert constant.dyncm2_to_MeVfm3 == type(float)
    assert constant.gcm3_to_MeVfm3 == type(float)
    assert constant.oneoverfm_MeV == type(float)

