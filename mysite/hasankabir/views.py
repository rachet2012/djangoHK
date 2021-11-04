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
        d_i_md  = float(form['d_i_m'].value()) 
        d_o_md = float(form['d_o_m'].value())
        hd = int(form['h'].value())
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
        
        data = Database()
        data.qu_liq_m3day = qu_liq_m3dayd
        data.d_i_m = d_i_md
        data.d_o_m = d_o_md
        data.h = hd
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

        ttt = schet(rp = rpd,qu_liq_r=qu_liq_m3dayd, wct_r=wctd, p_head_r = p_headd, t_head_r=(t_headd + 273),
                     d_i_r = d_i_md, d_o_r=d_o_md, absep_r = absepd, md1 = md1d,
                     md2 = md2d, md3 = md3d, tvd1 = tvd1d, tvd2 = tvd2d, tvd3=tvd3d)[0]
        data.results = ttt
        data.save()
        context [ 'result' ] = ttt

        lp = list(schet(rp = rpd,qu_liq_r=qu_liq_m3dayd, wct_r=wctd, p_head_r = p_headd, t_head_r=(t_headd + 273),
                     d_i_r = d_i_md, d_o_r=d_o_md, absep_r = absepd, md1 = md1d,
                     md2 = md2d, md3 = md3d, tvd1 = tvd1d, tvd2 = tvd2d, tvd3=tvd3d)[1])
        lp[0] =lp[0]*101325
        lp2 = []
        lh = [i for i in range(0, int(tvd3d+50), 300)]

        for i in lp:
            b = i/101325
            lp2.append(b)

        data = []
        data.append(Scattergl(y=lh, x=lp2, mode='lines',
                            line={'dash': 'solid', 'color': '#AF479D'}))
        layout = Layout(width=800, height=600, legend=dict(orientation="h", y=1.1),
            paper_bgcolor='rgba(0,0,0,0)', plot_bgcolor='rgba(0,0,0,0)',)
        figure = Figure(data=data, layout=layout)
        figure.update_xaxes(linewidth=2, linecolor='#A6A8AB', gridcolor='#A6A8AB')
        figure.update_yaxes(linewidth=2, linecolor='#A6A8AB', gridcolor='#A6A8AB')
        figure.update_yaxes(autorange="reversed")
        figure.update_xaxes(side = 'top')
        plot_fig = plot(figure, auto_open=False, output_type='div')
        context['plot'] = plot_fig
    else:
        form = DataForm()
        context['form'] = form

    return render(request, 'hasankabir/index.html', context)

# def index(request):
#     return render(request, 'hasankabir/index.html')


