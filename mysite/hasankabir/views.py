from django.shortcuts import render
from hasankabir.forms import DataForm
from hasankabir.models import Database
from hasankabir.djHK import HasanKabirAnn

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
        
        data = Database()
        data.qu_liq_m3day = qu_liq_m3dayd
        data.d_i_m = d_i_md
        data.d_o_m = d_o_md
        data.h = hd
        data.p_head = p_headd
        data.t_head = t_headd
        data.wct = wctd
        data.rp = rpd

        well = HasanKabirAnn(qu_liq_m3day=qu_liq_m3dayd, d_i_m=d_i_md,
                             d_o_m=d_o_md, h=hd, p_head=p_headd, 
                             t_head=t_headd, wct=wctd, rp= rpd)
        vr = float(well.down_pr())
        data.results = vr
        data.save()
        context [ 'result' ] = vr

    else:
        form = DataForm()
        context['form'] = form

    return render(request, 'hasankabir/index.html', context)

# def index(request):
#     return render(request, 'hasankabir/index.html')


