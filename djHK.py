import math as m
import numpy as np
from numpy.core.fromnumeric import put
import scipy.constants as CONST
import scipy.optimize as sp
from scipy.integrate import solve_ivp
from unifloc.pvt.fluid_flow import FluidFlow
import unifloc.pipe._friction as fr
from math import fabs
import pandas as pd
#Новый коэфициент трения, ошибка меньше. Старый убрать

class HasanKabirAnn(FluidFlow):
    """
    Класс для расчета градиента давления в затрубном пространстве
    Определяются структура потока, истинная концентрация газа, плотность смеси,
    градиенты на гравитацию, трение, ускорение.
    Вычисляется распеределение давление в затрубном пространстве.
    """
    def __init__(self, qu_gas_m3day:float = 10, qu_liq_m3day: float = 432, d_i_m: float = 73,
                 d_o_m: float = 142,theta: float = 90, h:float = 2400, p_head:float = 15, 
                 t_head:float = 20, wct:float = 0.1, abseps:float = 2.54, rp = 0) -> None:
        """
        :param qu_gas_m3day: дебит скважины по газу, м3/сут
        :param qu_liq_m3day: дебит скважины по жидкости, м3/сут
        :param rho_gas_kgm3: плотность газа, кг/м3
        :param d_i_m: внешний диаметр НКТ, мм
        :param d_o_m: внутренний диаметр ЭК, мм
        :param theta: угол наклона скважины
        :param h: глубина скважины, м 
        :param p_head: давление на устье скважины, атм
        :param t_head: температура на устье скважины, С
        :param wct: обводненность продукции, дол.ед
        :param abseps: абсолютная шероховатость стенок трубы, м
        :param rp: газовый фактор, м3/м3
        """
        self.qu_gas_m3sec = None #qu_gas_m3day / 86400
        self.qu_liq_m3sec = qu_liq_m3day / 86400
        self.rp = rp
        self.wct = wct

        self.p_head = p_head * (101325)
        self.t_head = t_head + 273

        self.abseps =abseps / 100000
        self.d_i_m = d_i_m / 1000
        self.d_o_m = d_o_m / 1000
        self.h = h
        self.theta = theta

        self.C0 = 1.2
        self.C1 = 1.15

        #рассчитанные
    
        self.flow_pattern_name = None

        self.rho_gas_kgm31 = None
        self.epsi = None
        self.rho_mix_kgm3 = None

        self.result_grad_pam = None

    def _calc_rash(self):
        """
        Метод для расчета общих параметров потока
        """
        f_m2 =  CONST.pi * ((self.d_o_m/2)**2 - (self.d_i_m/2)**2)
        self.d_equ_m = self.d_o_m - self.d_i_m
        self.vs_gas_msec = self.qu_gas_m3sec / f_m2
        self.vs_liq_msec = self.qu_liq_m3sec / f_m2
        self.v_mix_msec = (self.qu_gas_m3sec + self.qu_liq_m3sec) / f_m2
        self.k_ratio_d = self.d_i_m / self.d_o_m

    def calc_PVT(self, p, t):
        """
        Метод для расчета PVT-модели 
        :param p: текущее давление, Па 
        :param t: текущая температура, К
        """
        self.qu_oil = self.qu_liq_m3sec * (1-self.wct)
        pvt_model =  {"black_oil": {"gamma_gas": 0.7, "gamma_wat": 1, "gamma_oil": 0.8,
                                         "rp": self.rp,
                                         "oil_correlations":
                                          {"pb": "Standing", "rs": "Standing",
                                           "rho": "Standing","b": "Standing",
                                          "mu": "Beggs", "compr": "Vasquez"},
                            "gas_correlations": {"ppc": "Standing", "tpc": "Standing",
                                                  "z": "Dranchuk", "mu": "Lee"},
                             "water_correlations": {"b": "McCain", "compr": "Kriel",
                                                    "rho": "Standing", "mu": "McCain"}}}
        PVT = FluidFlow(self.qu_liq_m3sec, self.wct, pvt_model)
        PVT.calc_flow(p, t)
        self.qu_gas_m3sec = PVT.qg
        self.mu_gas_pasec = PVT.mug
        self.mu_liq_pasec = PVT.mul
        self.rho_gas_kgm31 = PVT.rg
        self.rho_liq_kgm3 = PVT.rl
        # self.rho_liq_kgm3 = PVT.rho_oil
        self.sigma_Nm = PVT.stlg
        self._calc_rash()

    def _mixture_velocity_Caetano(self, initial_f):
        """
        Метод для определения критической скорости dispersed bubbly flow
        TWO-PHASE FLOW IN VERTICAL AND INCLINED ANNULI
        """
        rho_m_rash_kgm3 = (self.qu_gas_m3sec / (self.qu_gas_m3sec + self.qu_liq_m3sec) 
                                * self.rho_gas_kgm31 + (1 - self.qu_gas_m3sec / (self.qu_gas_m3sec + 
                                self.qu_liq_m3sec)) * self.rho_liq_kgm3)
        friction_coeff = self._friction_coefv2(rho_m_rash_kgm3)
        right_part = 2 * (initial_f ** 1.2) * (friction_coeff** 0.4) * ((2 / self.d_equ_m) **
                     0.4) * ((self.rho_liq_kgm3 / self.sigma_Nm) ** 0.6 ) *(0.4 * self.sigma_Nm / (
                         (self.rho_liq_kgm3 - self.rho_gas_kgm31) * CONST.g) ** 0.5)

        left_part = 0.725 + 4.15 * (self.vs_gas_msec / initial_f) ** 0.5

        return right_part - left_part

    def _friction_coefv2(self, rho):
        """
        Метод для определения коэффициента трения
        :param rho: фактическая плотность ГЖС
        :return: коэффициент трения
        """
        self.frict = fr.Friction(self.d_o_m)
        eps = self.abseps / self.d_o_m
        mu_mix_pasec = (self.vs_liq_msec / self.v_mix_msec * self.mu_liq_pasec 
                            + self.vs_gas_msec / self.v_mix_msec * self.mu_gas_pasec)
        self.number_Re_s = self.frict.calc_n_re(rho, self.v_mix_msec, mu_mix_pasec)
        self.ff = self.frict.calc_norm_ff(self.number_Re_s, eps, 1)
        return self.ff

    def _calc_pattern(self):
        """
        Метод для расчета критических значений скоростей
        TWO-PHASE FLOW IN VERTICAL AND INCLINED ANNULI
        Upward Vertical Two-Phase Flow Through an Annulus—Part I for to dispersed
        """
        #bubble to slug transition [4]
        v_d_msec = 1.53 * (CONST.g * self.sigma_Nm * (self.rho_liq_kgm3 
                        - self.rho_gas_kgm31) / (self.rho_liq_kgm3)**2 ) ** 0.25
        self.vs_gas_bubble2slug_msec = ((self.C0 * self.vs_liq_msec + v_d_msec) / (4 
                                    - self.C0)) # * np.sin(self.theta * np.pi/180)

        #to annular transirion [17]
        self.vs_gas_2annular_msec = 3.1 * (self.sigma_Nm * CONST.g * (self.rho_liq_kgm3 - 
                                    self.rho_gas_kgm31) / (self.rho_gas_kgm31) ** 2) ** 0.25 + 1

        #bubble/slug to dispersed transition [6]
        self.v_m_krit2disp_msec = fabs(float(sp.fsolve(self._mixture_velocity_Caetano, 10, maxfev=10)))
        # self.v_m_krit2disp_msec = 5

        self._set_flow_pattrn()

    def _set_flow_pattrn(self):
        """
        Метод для определения структуры потока
        """
        if self.vs_gas_msec >= self.vs_gas_2annular_msec:
            self.flow_pattern = 4
            self.flow_pattern_name = 'Annular flow pattern - Кольцевой режим'
        elif self.vs_gas_msec <= self.vs_gas_bubble2slug_msec: #and self.v_mix_msec < self.v_m_krit2disp_msec:
            self.flow_pattern = 0
            self.flow_pattern_name = 'Bubble flow pattern - пузырьковый режим'
        elif self.vs_gas_msec >= self.vs_gas_bubble2slug_msec and (0.25 * self.vs_gas_msec) < 0.52 and self.v_mix_msec < self.v_m_krit2disp_msec:
            self.flow_pattern = 2
            self.flow_pattern_name = 'Slug flow pattern - Пробковый режим'
        elif self.vs_gas_msec >= self.vs_gas_bubble2slug_msec and (0.25 * self.vs_gas_msec) >= 0.52:
            self.flow_pattern = 3
            self.flow_pattern_name = 'Chug flow pattern - Вспененный режим'
        # elif self.vs_gas_msec <= self.vs_gas_bubble2slug_msec and self.v_mix_msec < self.v_m_krit2disp_msec:
        #     self.flow_pattern = 0
        #     self.flow_pattern_name = 'Bubble flow pattern - пузырьковый режим'
        elif self.v_mix_msec >= self.v_m_krit2disp_msec:
            self.flow_pattern = 1
            self.flow_pattern_name = 'Dispersed bubble flow pattern - дисперсионно-пузырьковый режим'       

    def _calc_bubbly(self) -> float:
        """
        Метод для расчета истинной объемной концентрации газа в bubbly flow
        """
        v_d_msec = (1.53 * (CONST.g * self.sigma_Nm * (self.rho_liq_kgm3 - self.rho_gas_kgm31)
                     / (self.rho_liq_kgm3)**2 ) ** 0.25) #3

        v_gas_msec = self.C0 * self.v_mix_msec + v_d_msec #1
        self.epsi = self.vs_gas_msec / v_gas_msec #2

    def _calc_slug_churn(self, C) -> float:
        """
        Метод для расчета истинной объемной концентрации газа в slug и churn flow
        :param С: параметр распределения газовой фазы в потоке
        """
        v_d_msec = (1.53 * (CONST.g * self.sigma_Nm * (self.rho_liq_kgm3 - self.rho_gas_kgm31) 
                        / (self.rho_liq_kgm3)**2 ) ** 0.25) #3

        self.v_dt_msec = 1.2 *(self.vs_gas_msec + self.vs_liq_msec) + 0.345 * (CONST.g * (self.d_i_m + self.d_o_m)) ** 0.5

        self.epsi_s =0.2 

        self.epsi_t = self.vs_gas_msec / (self.C1 * self.v_mix_msec + self.v_dt_msec) #7


        if self.vs_gas_msec > 0.4:

            self.len_s_m = 0.1 / self.epsi_s 
            self.epsi = (1 - self.len_s_m) * self.epsi_t + 0.1  #9a

        else:


            self.len_s_m = 0.25 * self.vs_gas_msec / self.epsi_s
            self.epsi = (1 - self.len_s_m) * self.epsi_t + 0.25 * self.vs_gas_msec #9b

    def _actual_film_length(self, initial_llf):
        """
        Метод для вычисления фактический длины пленки жидкости в slug/churn
        Уравнение  [44]
        
        """
        coef_b = (-2 * (1 - self.vs_gas_msec / self.v_dt_msec) * ((self.vs_gas_msec - self.v_gls #45
                    *(1 - self.h_ls) / self.v_dt_msec)) * self.len_ls + (2 / CONST.g) * (self.v_dt_msec 
                    - self.v_lls) ** 2 * self.h_ls **2) / (1 - self.vs_gas_msec / self.v_dt_msec) ** 2

        coef_c = (((self.vs_gas_msec - self.v_gls * (1 - self.h_ls)) / self.v_dt_msec * self.len_ls)
                        / (1 - self.vs_gas_msec / self.v_dt_msec)) ** 2
        return initial_llf ** 2 + coef_b * initial_llf + coef_c

    def _acceler_grad_p(self) :
        """
        Метод для нахождения градиента давления на ускорения в churn и slug flow
        [20,23,38,44,47] 
        :return: градиент давления на ускорение
        """
        self.h_ls = 1 - self.epsi_s
        h_lf = 1 - self.epsi_t
        # film_fick = 
        # h_lf = 
        self.len_ls = 16 * self.d_equ_m #38
        len_su = self.len_ls / self.len_s_m
        self.v_lls = ((self.vs_liq_msec + self.vs_gas_msec) - 1.53 * ((self.rho_liq_kgm3 - self.rho_gas_kgm31) #23
                    * CONST.g * self.sigma_Nm / (self.rho_liq_kgm3) ** 2) ** 0.25 * self.h_ls ** 0.5 * (1 - self.h_ls)) 
        self.v_gls = (1.53 * ((self.rho_liq_kgm3 - self.rho_gas_kgm31) * CONST.g * self.theta / (self.rho_liq_kgm3) #20
                        ** 2) ** 0.25 * self.h_ls ** 0.5) - self.v_lls
        try:
            self.act_len_lf = fabs(float(sp.broyden1(self._actual_film_length, 0.5)))
        except:
            self.act_len_lf = 0.001
        v_llf = (CONST.g * 2 * self.act_len_lf) ** 0.5 - self.v_dt_msec #47
        self.grad_p_acc = (self.rho_liq_kgm3 * (h_lf / len_su) * (v_llf - self.v_dt_msec) 
                    * (v_llf - self.v_lls))
        if self.grad_p_acc >= 0:
            var = self.grad_p_acc
        else:
            var = 0
        return var

    def _acceler_grad_p_annular(self):
        """
        Метод для расчета потерь на ускорения в кольцевом режиме потока
        """
        self.v_dt_msec = (0.345 + 0.1 * (self.d_i_m / self.d_o_m)) *((CONST.g * self.d_o_m * (
                        self.rho_liq_kgm3 - self.rho_gas_kgm31)/(self.rho_gas_kgm31)) ** 0.5)
        len_su= 1
        act_len_lf = len_su
        v_llf = (CONST.g * 2 * act_len_lf) ** 0.5 - self.v_dt_msec
        grad_p_acc_an = (self.rho_liq_kgm3 * (self.hl_total / len_su) * (v_llf - self.v_dt_msec) 
                    * v_llf)
        return grad_p_acc_an

    def _ratio_t(self):
        """
        Метод для вычисления отношения толщины пленок жидкости в кольцевом потоке [87,88]
        :return:
        
        """
        angle_wt_average = (1 / (1 - self.k_ratio_d ** 2) * (2 * m.asin(self.k_ratio_d) + 2 * #[88]
                             self.k_ratio_d  * (1 - self.k_ratio_d ** 2) ** 0.5 - CONST.pi * 
                            self.k_ratio_d ** 2))
        t_ratio = angle_wt_average / ((2 * CONST.pi - angle_wt_average) * self.k_ratio_d)
        return t_ratio

    def _calc_hl_total_annular(self):
        """
        Метод для вычисления концентрации жидкости в потоке кольцевой структуры[77]
        Допустил, что толщина пленки жидкости на внешней(или внутренней) трубе известна()
        """
        delta_o = 0.005
        delta_i = delta_o * self._ratio_t()
        phi = (10 ** 4 * self.vs_gas_msec * self.mu_gas_pasec / self.sigma_Nm * (self.rho_gas_kgm31 #[79]
                    / self.rho_liq_kgm3) ** 0.5)
        fe = 1 - m.exp((phi - 1.5) * (-0.125))
        self.hl_total = (4 / ((self.d_o_m) * (1 - self.k_ratio_d ** 2)) * (delta_o * (1 - delta_o / self.d_o_m) #[77]
                        + delta_i * self.k_ratio_d * (1 + delta_i / self.d_i_m) + self.vs_liq_msec * fe 
                        / ((self.vs_liq_msec * fe + self.vs_gas_msec) * (1 - self.k_ratio_d ** 2)) * 
                        (1 - self.k_ratio_d ** 2 - 4 * delta_o / self.d_o_m * (1 - delta_o / self.d_o_m) 
                        - 4 * delta_i *self.k_ratio_d / self.d_o_m * (1 + delta_i / self.d_i_m))))
        self.epsi = 1 - self.hl_total 
          
    def calc_rho_mix(self):
        """
        Метод для расчета плотности смеси
        """
        self._calc_pattern()
        if self.flow_pattern == 0 or self.flow_pattern == 1:
            self._calc_bubbly()
        elif self.flow_pattern == 2:
            self._calc_slug_churn(self.C0)
        elif self.flow_pattern == 3: 
            C1 = 1.15
            self._calc_slug_churn(C1)
        elif self.flow_pattern == 4:
            self._calc_hl_total_annular()
        
        self.rho_mix_kgm3 = self.rho_liq_kgm3 * (1 - self.epsi) + self.rho_gas_kgm31 * self.epsi

    def fanning_f(self,rho) :
        self.frict = fr.Friction(self.d_o_m)
        mu_mix_pasec = (self.vs_liq_msec / self.v_mix_msec * self.mu_liq_pasec 
                            + self.vs_gas_msec / self.v_mix_msec * self.mu_gas_pasec)
        Re = self.frict.calc_n_re(rho, self.v_mix_msec, mu_mix_pasec)
        abs_rough = self.abseps * 3.281
        d = self.d_equ_m * 3.281
        f_lam = 16 / Re
        inv_sqrt_f = -4 * np.log10(0.2698 * (abs_rough / d) - (5.0452 / Re) * np.log10(.3539 * (abs_rough / d)**1.1098 + 5.8506 / Re**.8981))
        f_turb = (1 / inv_sqrt_f)**2
    
        if(Re <= 2000):
            self.fri = f_lam
        elif((Re > 2000) & (Re < 4000)):
            self.fri = (f_lam * (4000 - Re) + f_turb * (Re - 2000)) / 2000
        elif(Re >= 4000):
            self.fri = f_turb

        return self.fri

    def calc_pressure_gradient(self, p, t):
        """
        Метод для расчета градиента давления
        Upward Vertical Two-Phase Flow Through an Annulus—Part II
        :param p: текущее давление, Па 
        :param t: текущая температура, К
        :return: суммарный градиент давления
        
        """
        self.calc_PVT(p, t)
        self.calc_rho_mix()

        if self.flow_pattern == 0: #[5-14]
            self.density_grad_pam = self.rho_mix_kgm3 * CONST.g * np.sin(self.theta * np.pi/180)

            # friction_coeff_s = self._friction_coefv2(self.rho_mix_kgm3)
            friction_coeff_s = self.fanning_f(self.rho_mix_kgm3)
            self.friction_grad_pam = ((4 * friction_coeff_s / (self.d_equ_m ) 
                                     * self.v_mix_msec ** 2 / 2) * self.rho_mix_kgm3)

            self.acceleration_grad_pam = 0

        elif self.flow_pattern == 1: #[15-16]
            self.density_grad_pam = self.rho_mix_kgm3 * CONST.g * np.sin(self.theta * np.pi/180)

            friction_coeff_s = self._friction_coefv2(self.rho_mix_kgm3)
            self.friction_grad_pam = (4 * friction_coeff_s / (self.d_equ_m ) 
                                     * self.v_mix_msec ** 2 / 2) * self.rho_mix_kgm3

            self.acceleration_grad_pam = 0
            
        elif self.flow_pattern == 2 or self.flow_pattern == 3: #предположил что для slug и churn одна методика. Концентрацию воды нашел как 1 - epsi
            self.rho_slug_kgm3 = self.rho_gas_kgm31 * self.epsi_s + self.rho_liq_kgm3 * (1 - self.epsi_s) # В соответствии c [51] 

            self.density_grad_pam = self.rho_slug_kgm3 * CONST.g  * self.len_s_m #[50]
            
            friction_coeff_s = self.fanning_f(self.rho_slug_kgm3)
            self.friction_grad_pam = ((2 * friction_coeff_s / self.d_equ_m * self.rho_slug_kgm3) #[53]
                                     * (self.vs_gas_msec + self.vs_liq_msec) **2 * self.len_s_m)

            self.acceleration_grad_pam = self._acceler_grad_p() 

        elif self.flow_pattern == 4:# над ускорением подумать
            self.density_grad_pam = self.rho_mix_kgm3 * CONST.g * np.sin(self.theta * np.pi/180)

            friction_coeff_s = self._friction_coefv2(self.rho_mix_kgm3)
            self.friction_grad_pam = (4 * friction_coeff_s / (self.d_o_m - self.d_i_m) 
                                     * self.v_mix_msec ** 2 / 2) * self.rho_mix_kgm3

            self.acceleration_grad_pam = self._acceler_grad_p_annular()

        self.result_grad_pam = self.friction_grad_pam  + self.density_grad_pam + self.acceleration_grad_pam

        print(self.flow_pattern_name)

        return self.result_grad_pam

    def _grad_func(self, h, pt):
        """
        Функция для интегрирования 
        :param pt: давление и температура, Па,К
        :h: переменная интегрирования
        """ 
        dp_dl = self.calc_pressure_gradient(pt[0], pt[1]) 
        dt_dl = 0.03
        return dp_dl, dt_dl

    def func_p_list(self):
        """
        Метод для интегрирования градиента давления(расчет сверху вниз)
        :return: распределение давления и температуры по стволу скважины
        """
        p0,t0 = self.p_head, self.t_head
        h0 = 0
        h1 = self.h
        steps = [i for i in range(h0, h1+50, 300)]
        sol = solve_ivp(self._grad_func, 
            t_span=(h0, h1), 
            y0=[p0, t0], 
            t_eval=steps,
            max_step=2050,
            rtol=1e-6,
            # event = enent
        ) 
        return sol.y, 



if __name__ == '__main__':

    #ТЕСТ
    p1 = []
    rbb1 =[]
    p2 = []
    rbb2 =[]
    p3 = []
    rbb3 =[]
    p4 = []
    rbb4 =[]
    
    for i in range(100,110, 10):
        rb =i
        test2 = HasanKabirAnn(rp =rb, qu_liq_m3day=285,wct = 0.5)
        vr = test2.func_p_list()
        vr1 = vr[0]
        vr2 = vr1[0]
        vr3= vr2[-1]
        rbb1.append(rb)
        p1.append(vr3)
        print(vr3/101325)
        
        # print(vr2)
    # for i in range(0,200, 10):
    #     rb =i
    #     test2 = HasanKabirAnn(rp =rb, qu_liq_m3day=285,wct=0)
    #     vr = test2.func_p_list()
    #     vr1 = vr[0]
    #     vr2 = vr1[0]
    #     vr3= vr2[-1] / 101325
    #     rbb1.append(rb)
    #     p1.append(vr3)
    # for i in range(0,200, 10):
    #     rb =i
    #     test2 = HasanKabirAnn(rp =rb, qu_liq_m3day=285,wct=0.1)
    #     vr = test2.func_p_list()
    #     vr1 = vr[0]
    #     vr2 = vr1[0]
    #     vr3= vr2[-1] / 101325
    #     rbb2.append(rb)
    #     p2.append(vr3)
    # for i in range(0,200, 10):
    #     rb =i
    #     test2 = HasanKabirAnn(rp =rb, qu_liq_m3day=285,wct=0.25)
    #     vr = test2.func_p_list()
    #     vr1 = vr[0]
    #     vr2 = vr1[0]
    #     vr3= vr2[-1] / 101325
    #     rbb3.append(rb)
    #     p3.append(vr3)
    # for i in range(0,200, 10):
    #     rb =i
    #     test2 = HasanKabirAnn(rp =rb, qu_liq_m3day=285,wct=0.4)
    #     vr = test2.func_p_list()
    #     vr1 = vr[0]
    #     vr2 = vr1[0]
    #     vr3= vr2[-1] / 101325
    #     rbb4.append(rb)
    #     p4.append(vr3)

    # df1 = pd.DataFrame({'GOR': rbb1,
    #                'p down': p1})

    # df2 = pd.DataFrame({'GOR': rbb2,
    #                'p down': p2})
                   
    # df3 = pd.DataFrame({'GOR': rbb3,
    #                'p down': p3})

    # df4 = pd.DataFrame({'GOR': rbb4,
    #                'p down': p4})    

    # salary_sheets = {'q285,wct0': df1, 'q285,wct0.1': df2, 'q285,wct0.25':df3, 'q285,wct0.4':df4}
    # writer = pd.ExcelWriter('./test4.xlsx', engine='xlsxwriter')

    # for sheet_name in salary_sheets.keys():
    #     salary_sheets[sheet_name].to_excel(writer, sheet_name=sheet_name, index=False)

    # writer.save()

#     df = pd.DataFrame({'GOR': rbb,
#                    'p down': p})
#     df.to_excel('./test.xlsx')