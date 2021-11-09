from django import forms

class DataForm(forms.Form):
    qu_liq_m3day = forms.FloatField(label='Дебит скважины по жидкости, м3/сут:', initial=300)
    wct = forms.FloatField(label='Обводненность продукции, дол.ед:', initial=0.6)
    rp = forms.FloatField(label='Газовый фактор, м3/м3:', initial=100)
    p_head = forms.FloatField(label='Давление на устье скважины, атм:', initial=15)
    t_head = forms.FloatField(label='Температура на устье скважины, С:', initial=20)
    gamma_gas = forms.FloatField(label='Относительная плотность газа, дол.ед:', initial=0.7)
    gamma_wat = forms.FloatField(label='Относительная плотность воды, дол.ед:', initial=1)
    gamma_oil = forms.FloatField(label='Относительная плотность нефти, дол.ед:', initial=0.8)
    pb = forms.FloatField(label='Давление насыщения, атм:', initial=98)
    t_res = forms.FloatField(label='Пластовая температура, C:', initial=92)
    rsb = forms.FloatField(label='Газосодержание нефти при давлении насыщения, м3/м3:', initial=50)
    muob = forms.FloatField(label='Вязкость нефти при давлении насыщения, сПз:', initial=0.5)
    bob = forms.FloatField(label='Объемный коэффициент нефти при давлении насыщения, м3/м3:', initial=1.5)
    absep = forms.FloatField(label='Абсолютная шероховатость стеном трубы, м*10^-5 :', initial=2.54)
    md1 = forms.IntegerField(label='MD p1, м :', initial=1400)
    md2 = forms.IntegerField(label='MD p2, м :', initial=1800)
    md3 = forms.IntegerField(label='MD p3, м :', initial=2400)
    tvd1 = forms.IntegerField(label='TVD p1, м :', initial=1400)
    tvd2 = forms.IntegerField(label='TVD p2, м :', initial=1800)
    tvd3 = forms.IntegerField(label='TVD p3, м :', initial=2400)
    d_o_1 = forms.FloatField(label='Внутренний диаметр ЭК в MD1, мм:', initial=142)
    d_o_2 = forms.FloatField(label='Внутренний диаметр ЭК в MD2, мм:', initial=142)
    d_o_3 = forms.FloatField(label='Внутренний диаметр ЭК в MD3, мм:', initial=142)
    d_i_1 = forms.FloatField(label='Внешний диаметр НКТ  в MD1, мм:', initial=73)
    d_i_2 = forms.FloatField(label='Внешний диаметр НКТ  в MD2, мм:', initial=73)
    d_i_3 = forms.FloatField(label='Внешний диаметр НКТ  в MD3, мм:', initial=73)






