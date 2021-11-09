from collections import namedtuple
from django.shortcuts import render
from hasankabir.forms import DataForm
from hasankabir.models import Database
from hasankabir.djHK import *
from plotly.offline import plot
from plotly.graph_objects import Scattergl, Layout, Figure

def index(request):
    context = {}
    if request.method == 'POST':
        form = DataForm(request.POST)
        context['form'] = form

        qu_liq_m3dayd = float(form['qu_liq_m3day'].value())
        p_headd = float(form['p_head'].value())
        t_headd = float(form['t_head'].value())
        wctd = float(form['wct'].value())
        rpd = float(form['rp'].value())
        absepd = float(form['absep'].value())
        md1d = float(form['md1'].value())
        md2d = float(form['md2'].value())
        md3d = float(form['md3'].value())
        tvd1d = float(form['tvd1'].value())
        tvd2d = float(form['tvd2'].value())
        tvd3d = float(form['tvd3'].value())
        gamma_gas = float(form['gamma_gas'].value())
        gamma_wat = float(form['gamma_wat'].value())
        gamma_oil = float(form['gamma_oil'].value())
        pb = float(form['pb'].value())
        t_res = float(form['t_res'].value())
        rsb = float(form['rsb'].value())
        muob = float(form['muob'].value())
        bob = float(form['bob'].value())
        d_o_1 = float(form['d_o_1'].value())
        d_o_2 = float(form['d_o_2'].value())
        d_o_3 = float(form['d_o_3'].value())
        d_i_1 = float(form['d_i_1'].value())
        d_i_2 = float(form['d_i_2'].value())
        d_i_3 = float(form['d_i_3'].value())


        data = Database()
        data.qu_liq_m3day = qu_liq_m3dayd
        data.p_head = p_headd
        data.t_head = t_headd
        data.wct = wctd
        data.rp = rpd
        data.absep = absepd
        data.md1 = md1d
        data.md2 = md2d
        data.md3 = md3d
        data.tvd1 = tvd1d
        data.tvd2 = tvd2d
        data.tvd3 = tvd3d
        data.gamma_gas = gamma_gas
        data.gamma_wat = gamma_wat
        data.gamma_oil = gamma_oil
        data.pb = pb
        data.t_res = t_res
        data.rsb = rsb
        data.muob = muob
        data.bob = bob
        data.d_o_1 = d_o_1
        data.d_o_2 = d_o_2
        data.d_o_3 = d_o_3
        data.d_i_1 = d_i_1
        data.d_i_2 = d_i_2
        data.d_i_3 = d_i_3


        ttt = schet(rp = rpd,qu_liq_r=qu_liq_m3dayd, wct_r=wctd, p_head_r = (p_headd*101325), t_head_r=(t_headd + 273),
                    absep_r = absepd, md1 = md1d,
                     md2 = md2d, md3 = md3d, tvd1 = tvd1d, tvd2 = tvd2d, tvd3=tvd3d,
                     gamma_gas=gamma_gas, gamma_oil = gamma_oil,gamma_wat= gamma_wat, pb = (pb*101325),t_res=(t_res+273),rsb=rsb,
                     muob=muob,bob=bob, d_o_1 = d_o_1,d_o_2 = d_o_2,d_o_3 = d_o_3,d_i_1=d_i_1,d_i_2=d_i_2,d_i_3=d_i_3)[0]
        data.results = ttt
        data.save()
        context [ 'result' ] = ttt

        lp = list(schet(rp = rpd,qu_liq_r=qu_liq_m3dayd, wct_r=wctd, p_head_r = (p_headd*101325), t_head_r=(t_headd + 273),
                     absep_r = absepd, md1 = md1d,
                     md2 = md2d, md3 = md3d, tvd1 = tvd1d, tvd2 = tvd2d, tvd3=tvd3d,
                     gamma_gas=gamma_gas, gamma_oil = gamma_oil,gamma_wat= gamma_wat, pb = (pb*101325),t_res=(t_res+273),rsb=rsb,
                     muob=muob,bob=bob,d_o_1 = d_o_1,d_o_2 = d_o_2,d_o_3 = d_o_3,d_i_1=d_i_1,d_i_2=d_i_2,d_i_3=d_i_3)[1])
        lp2 = []
        lh = [i for i in range(0, int(tvd3d+50), 50)]

        for i in lp:
            b = i/101325
            lp2.append(b)

        data = []
        data.append(Scattergl(y=lh, x=lp2, mode='lines + markers',
                            line={'dash': 'solid', 'color': '#4B0082'},
                            marker = {'color': '#FF0000',} ))
        layout = Layout(width=800, height=600, legend=dict(orientation="h", y=0),
            paper_bgcolor='rgba(0,0,0,0)', plot_bgcolor='rgba(0,0,0,0)',)
        figure = Figure(data=data, layout=layout)
        figure.update_xaxes(linewidth=2, linecolor='#A6A8AB', gridcolor='#A6A8AB', range=[0, ttt+50])
        figure.update_yaxes(linewidth=2, linecolor='#A6A8AB', gridcolor='#A6A8AB', range=[0, tvd3d+50])
        figure.update_yaxes(autorange="reversed", zeroline=True, zerolinewidth=2, zerolinecolor='LightPink')
        figure.update_xaxes(side = 'top', zeroline=True, zerolinewidth=2, zerolinecolor='LightPink')
        figure.update_layout(xaxis_title="Давление, атм",
                  yaxis_title="Глубина, м", title ='КРД')
        figure.update_traces(hoverinfo="all", hovertemplate="Давление, атм: %{x}<br>Глубина,м: %{y}")
        plot_fig = plot(figure, auto_open=False, output_type='div')
        context['plot'] = plot_fig
    else:
        form = DataForm()
        context['form'] = form

    return render(request, 'hasankabir/index.html', context)

# def index(request):
#     return render(request, 'hasankabir/index.html')


