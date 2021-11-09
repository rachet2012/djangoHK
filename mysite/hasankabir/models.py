from django.db import models


class Database(models.Model):
    """
    :param qu_liq_m3day: дебит скважины по жидкости, м3/сут
    :param d_i_m: внешний диаметр НКТ, мм
    :param d_o_m: внутренний диаметр ЭК, мм
    :param h: глубина скважины, м
    :param p_head: давление на устье скважины, атм
    :param t_head: температура на устье скважины, С
    :param wct: обводненность продукции, дол.ед
    """

    qu_liq_m3day = models.FloatField()
    p_head = models.FloatField()
    t_head = models.FloatField()
    wct = models.FloatField()
    rp = models.FloatField()
    absep = models.FloatField()
    md1 = models.FloatField()
    md2 = models.FloatField()
    md3 = models.FloatField()
    tvd1 = models.FloatField()
    tvd2 = models.FloatField()
    tvd3 = models.FloatField()
    gamma_gas = models.FloatField()
    gamma_wat = models.FloatField()
    gamma_oil = models.FloatField()
    pb = models.FloatField()
    t_res = models.FloatField()
    rsb = models.FloatField()
    muob = models.FloatField()
    bob = models.FloatField()
    d_o_1 = models.FloatField()
    d_o_2 = models.FloatField()
    d_o_3 = models.FloatField()
    d_i_1 = models.FloatField()
    d_i_2 = models.FloatField()
    d_i_3 = models.FloatField()
    results = models.CharField(max_length=1000000)


    def __str__(self):
        return self.results

