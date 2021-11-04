from django import forms

class DataForm(forms.Form):
    qu_liq_m3day = forms.FloatField(label='Дебит скважины по жидкости, м3/сут:', initial=300)
    wct = forms.FloatField(label='Обводненность продукции, дол.ед:', initial=0.6)
    rp = forms.FloatField(label='Газовый фактор, м3/м3:', initial=100)
    p_head = forms.FloatField(label='Давление на устье скважины, атм:', initial=15)
    t_head = forms.FloatField(label='Температура на устье скважины, С:', initial=20)
    d_i_m = forms.FloatField(label='Внешний диаметр НКТ, мм:', initial=73)
    d_o_m = forms.FloatField(label='Внутренний диаметр ЭК, мм:', initial=142)
    absep = forms.IntegerField(label='Абсолютная шероховатость стеном трубы, м*10^-5 :', initial=2.54)
    md1 = forms.IntegerField(label='MD p1, м :', initial=1400)
    md2 = forms.IntegerField(label='MD p2, м :', initial=1800)
    md3 = forms.IntegerField(label='MD p3, м :', initial=2400)
    tvd1 = forms.IntegerField(label='TVD p1, м :', initial=1400)
    tvd2 = forms.IntegerField(label='TVD p2, м :', initial=1800)
    tvd3 = forms.IntegerField(label='TVD p3, м :', initial=2400)    





    
