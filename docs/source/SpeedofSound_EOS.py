import numpy as np
import math
from scipy import optimize
from TOVsolver.unit import g_cm_3, dyn_cm_2, e0
from scipy.integrate import cumulative_simpson


#####################  SPEED-OF-SOUND-EOS  BEGIN  ##############################


class compute_EOS:
    def __init__(self, x_last, y_last, dydx_last, enablePTcheck=False) -> None:
        """
        Speed of sound EOS class.
        See https://arxiv.org/abs/1812.08188 for details.
        Since the speed of sound core eos need to fit the last point and derivative of outer eos,
        the initial values are provided here.

        Args:
            x_last (float): Last energy density value.
            y_last (float): Last pressure value.
            dydx_last (float): Last derivative value.
            enablePTcheck (bool, optional): Enable phase transition check. Defaults to False.
        """
        self.x_last = x_last
        self.y_last = y_last
        self.dydx_last = dydx_last
        self.enablePTcheck = enablePTcheck

    def cs2(self, x, a):
        a1, a2, a3, a4, a5, a6 = a
        ret = (
            a1 * np.exp(-0.5 * (((x - a2) / a3) ** 2))
            + a6
            + (1 / 3 - a6) / (1 + np.exp(-a5 * (x - a4)))
        )
        return np.clip(ret, 0, 1)  # requirement 2 and 4

    def cal_a6(self, a1, a2, a3, a4, a5):
        A = a1 * np.exp(-0.5 * (((self.x_last - a2) / a3) ** 2)) + (1 / 3) / (
            1 + np.exp(-a5 * (self.x_last - a4))
        )
        B = 1 - 1 / (1 + np.exp(-a5 * (self.x_last - a4)))
        # cs2 = A + B * a6 = dydx
        return (self.dydx_last - A) / B

    def uniform_prior(self, a, b, r):
        return a + (b - a) * r

    def gen_a(self, cube5_):
        """
        Generate the EOS parameters from the given cube.

        Args:
            cube5_ (array): The input cube. The shape is (5,). The values are in the range of [0, 1].

        Returns:
            tuple: The EOS parameters.
        """
        a1 = self.uniform_prior(0.1, 1.5, cube5_[0])
        a2 = self.uniform_prior(1.5, 12, cube5_[1])
        a3 = self.uniform_prior(0.05 * a2, 2 * a2, cube5_[2])
        a4 = self.uniform_prior(1.5, 37, cube5_[3])
        a5 = self.uniform_prior(0.1, 1, cube5_[4])
        a6 = self.cal_a6(a1, a2, a3, a4, a5)
        return (a1, a2, a3, a4, a5, a6)

    def check_a(self, a):
        a1, a2, a3, a4, a5, a6 = a

        # PT requirement
        if self.enablePTcheck and a6 > 0:
            return False

        # requirement 5
        if np.any(self.cs2(np.linspace(0.5 * e0, 1.5 * e0, 16), a) > 0.163):
            return False

        # requirement 3
        if a6 > 1 / 3:
            return False
        if np.any(self.cs2(np.linspace(49 * e0, 50 * e0, 16), a) > 1 / 3):
            return False

        return True

    def cal_core_p(self, core_e, a):
        """
        Calculate the core pressure.

        Args:
            core_e (array): The core energy density.
            a (array): The EOS parameters.

        Returns:
            array: The core pressure.
        """
        core_p = cumulative_simpson(self.cs2(core_e, a), x=core_e, initial=self.y_last)
        return core_p
        # fix check dpde < 0 issue
        # dif_p = np.diff(core_p, prepend=self.y_last)
        # dif_p = np.maximum(dif_p, 0)
        # return np.cumsum(dif_p)


#####################  SPEED-OF-SOUND-EOS  END  ###########################  ###
