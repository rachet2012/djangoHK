from django import forms

class DataForm(forms.Form):
    qu_liq_m3day = forms.FloatField(label='Дебит скважины по жидкости, м3/сут:', initial=285)
    d_i_m = forms.FloatField(label='Внешний диаметр НКТ, мм:', initial=73)
    d_o_m = forms.FloatField(label='Внутренний диаметр ЭК, мм:', initial=142)
    h = forms.IntegerField(label='Глубина скважины, м :', initial=2400)
    p_head = forms.FloatField(label='Давление на устье скважины, атм:', initial=15)
    t_head = forms.FloatField(label='Температура на устье скважины, С:', initial=20)
    wct = forms.FloatField(label='Обводненность продукции, дол.ед:', initial=0.5)
    rp = forms.FloatField(label='Газовый фактор, м3/м3:', initial=100)



    
