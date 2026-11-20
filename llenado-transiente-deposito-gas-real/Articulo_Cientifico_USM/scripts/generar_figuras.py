"""Generación reproducible de figuras del manuscrito ICF130 (solo CSV del proyecto)."""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
D=Path(__file__).resolve().parents[1]/'datos'
F=Path(__file__).resolve().parents[1]/'figuras'
F.mkdir(exist_ok=True)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9.1,'axes.titlesize':10.2,'axes.labelsize':9.1,'axes.linewidth':0.8,'axes.spines.top':False,'axes.spines.right':False,'xtick.labelsize':8.2,'ytick.labelsize':8.2,'legend.fontsize':7.6,'figure.dpi':160,'savefig.bbox':'tight','pdf.fonttype':42})
COL={'Ideal':'#215a95','Van der Waals':'#ba6334','Peng–Robinson':'#0b8379'}
STYL={'Ideal':'-','Van der Waals':'--','Peng–Robinson':'-.'}
models=list(COL)
read=lambda name:pd.read_csv(D/name,encoding='utf-8-sig')
R={f'{v}_{m}':read(f'{v}_{m}.csv') for v in ['V1','V2','V3'] for m in ['ideal','vdw','pr']}
B=read('casos_base.csv'); C=read('comparacion_controlada_V1_V2.csv'); M=read('sensibilidad_masa_pared.csv'); S=read('sensibilidad_pared_3EOS.csv'); H=read('propiedades_contra_HEOS.csv'); GRID=read('V3_malla_suministro.csv');ISO=read('V2_limite_isotermo.csv'); RK=read('orden_RK4_paso_fijo_V3.csv')
keys=dict(zip(models,['ideal','vdw','pr']))
def styling(ax):
 ax.grid(color='#e5e9ee',lw=.55,alpha=.9)
 ax.set_axisbelow(True)
def finish(fig,name):
 fig.tight_layout(pad=1.0)
 fig.savefig(F/(name+'.pdf'))
 fig.savefig(F/(name+'.png'),dpi=210)
 plt.close(fig)

# F01 V1 comparison
fig,axs=plt.subplots(1,2,figsize=(7.25,2.75))
for m in models:
 df=R['V1_'+keys[m]]; axs[0].plot(df['t_s'],df['P_bar'],c=COL[m],label=m,lw=1.75);axs[1].plot(df['t_s'],df['T_K'],c=COL[m],lw=1.75)
axs[0].set(xlabel='Tiempo [s]',ylabel='Presión [bar]',title='(a) Presión: V1, caudal fijo');axs[1].set(xlabel='Tiempo [s]',ylabel='Temperatura del gas [K]',title='(b) Temperatura: V1')
for ax in axs: styling(ax)
fig.legend(*axs[0].get_legend_handles_labels(),loc='lower center',bbox_to_anchor=(.5,-.09),ncol=3,frameon=False)
fig.subplots_adjust(bottom=.17)
finish(fig,'01_v1_trayectorias')
# F02 flow
fig,axs=plt.subplots(1,2,figsize=(7.25,2.75))
for m in models:
 df=R['V2_'+keys[m]];axs[0].plot(df['P_bar'],df['mdot_g_s'],label=m,c=COL[m],lw=1.8);axs[1].plot(df['t_s'],df['mdot_g_s'],c=COL[m],lw=1.8)
axs[0].axvline(105.66,c='0.3',ls=':',lw=1.1);axs[0].text(105.66,11,'105,66 bar',rotation=90,ha='right',fontsize=7,color='0.3')
axs[0].set(xlabel='Presión interna [bar]',ylabel='Caudal [g/s]',title='(a) Ley común en función de P');axs[1].set(xlabel='Tiempo [s]',ylabel='Caudal [g/s]',title='(b) Caudal a lo largo del llenado')
for ax in axs: styling(ax)
fig.legend(*axs[0].get_legend_handles_labels(),loc='lower center',bbox_to_anchor=(.5,-.09),ncol=3,frameon=False)
finish(fig,'02_v2_caudal')
# F03 V3 wall and gas
fig,axs=plt.subplots(1,2,figsize=(7.25,2.75))
for m in models:
 df=R['V3_'+keys[m]];axs[0].plot(df.t_s,df.T_K,c=COL[m],label=m,lw=1.8);axs[1].plot(df.t_s,df.Tw_K,c=COL[m],label=m,lw=1.8)
axs[0].set(xlabel='Tiempo [s]',ylabel='T gas [K]',title='(a) Temperatura del gas');axs[1].set(xlabel='Tiempo [s]',ylabel='T pared [K]',title='(b) Temperatura de pared')
for ax in axs: styling(ax)
fig.legend(*axs[0].get_legend_handles_labels(),loc='lower center',bbox_to_anchor=(.5,-.09),ncol=3,frameon=False)
finish(fig,'03_v3_temperaturas')
# F04 V3 pressure mass Z
fig,axs=plt.subplots(1,3,figsize=(7.4,2.75))
for m in models:
 df=R['V3_'+keys[m]]
 for ax,col in zip(axs,['P_bar','masa_kg','Z']):ax.plot(df.t_s,df[col],c=COL[m],label=m,lw=1.75)
for ax,y,title in zip(axs,['P [bar]','m [kg]','Z [–]'],['(a) Presión','(b) Masa','(c) Compresibilidad']):ax.set(xlabel='Tiempo [s]',ylabel=y,title=title);styling(ax)
fig.legend(*axs[0].get_legend_handles_labels(),loc='lower center',bbox_to_anchor=(.5,-.08),ncol=3,frameon=False)
finish(fig,'04_v3_P_m_Z')
# F05 controlled V1 V2
fig,axs=plt.subplots(1,2,figsize=(7.25,2.75));xx=np.arange(3);w=.25
for i,m in enumerate(models):
 df=C[C.eos==m].sort_values('Pobj_bar')
 axs[0].plot(df.Pobj_bar,df.diferencia_s,marker='o',ms=5,label=m,lw=1.7,c=COL[m]);
 axs[1].bar(xx[i]-w/2,df.iloc[-1].t_V1_s,width=w,color=COL[m],alpha=.53)
 axs[1].bar(xx[i]+w/2,df.iloc[-1].t_V2_s,width=w,color=COL[m],alpha=1)
axs[0].set(xlabel='Presión objetivo [bar]',ylabel='Aumento V2 − V1 [s]',title='(a) Contraste controlado');axs[0].set_xticks([100,190]);axs[1].set(xlabel='EOS (V1 claro / V2 oscuro)',ylabel='Tiempo hasta 190 bar [s]',title='(b) Mismo objetivo, distinto caudal');axs[1].set_xticks(xx,labels=['Ideal','VdW','PR']);
for ax in axs:styling(ax)
fig.legend(*axs[0].get_legend_handles_labels(),loc='lower center',bbox_to_anchor=(.5,-.085),ncol=3,frameon=False)
finish(fig,'05_control_V1_V2')
# F06 thermal sensitivity UA
fig,axs=plt.subplots(1,2,figsize=(7.25,2.75)); sub=S[S.parametro=='UA_gas_pared']
for m in models:
 d=sub[sub.modelo==m].sort_values('valor');axs[0].plot(d.valor,d.tf_s,marker='o',ms=4,c=COL[m],lw=1.6,label=m);axs[1].plot(d.valor,d.Tmax_K,marker='o',ms=4,c=COL[m],lw=1.6)
axs[0].set(xlabel=r'$UA_{gw}$ [W/K]',ylabel='Tiempo [s]',title='(a) Duración hasta 190 bar');axs[1].set(xlabel=r'$UA_{gw}$ [W/K]',ylabel='Temperatura máxima [K]',title='(b) Máximo térmico')
for ax in axs:styling(ax)
fig.legend(*axs[0].get_legend_handles_labels(),loc='lower center',bbox_to_anchor=(.5,-.085),ncol=3,frameon=False)
finish(fig,'06_sensibilidad_UA')
# F07 wall mass
fig,axs=plt.subplots(1,2,figsize=(7.25,2.75))
for m in models:
 d=M[M.eos==m].sort_values('masa_efectiva_kg'); axs[0].plot(d.masa_efectiva_kg,d.temperatura_final_K,c=COL[m],marker='o',label=m,lw=1.6);axs[1].plot(d.masa_efectiva_kg,d.t_final_s,c=COL[m],marker='o',lw=1.6)
axs[0].set(xlabel='Masa térmica efectiva [kg]',ylabel='Temperatura final [K]',title='(a) Sensibilidad de temperatura');axs[1].set(xlabel='Masa térmica efectiva [kg]',ylabel='Tiempo final [s]',title='(b) Sensibilidad de duración')
for ax in axs:styling(ax)
fig.legend(*axs[0].get_legend_handles_labels(),loc='lower center',bbox_to_anchor=(.5,-.085),ncol=3,frameon=False)
finish(fig,'07_masa_pared')
# F08 heatmap supply grid
fig,axs=plt.subplots(1,2,figsize=(7.25,2.75))
for ax,value,name in [(axs[0],'tf_s','Tiempo PR − ideal [s]'),(axs[1],'Tf_K','Temperatura PR − ideal [K]')]:
 pv=GRID.pivot_table(index=['Ps_bar','Ts_K'],columns='modelo',values=value)
 d=(pv['Peng–Robinson']-pv['Ideal']).unstack('Ts_K').sort_index().sort_index(axis=1)
 import matplotlib.colors as colors
 if value=='tf_s': cmap='RdBu_r'; norm=colors.TwoSlopeNorm(vcenter=0,vmin=-16,vmax=9)
 else: cmap='PuBu';norm=None
 im=ax.imshow(d.values,origin='lower',aspect='auto',cmap=cmap,norm=norm)
 ax.set_xticks(range(len(d.columns)),[str(int(x)) for x in d.columns]);ax.set_yticks(range(len(d.index)),[str(int(x)) for x in d.index]);ax.set(xlabel=r'$T_s$ [K]',ylabel=r'$P_s$ [bar]',title=name)
 for (i,j),v in np.ndenumerate(d.values):ax.text(j,i,f'{v:+.1f}',ha='center',va='center',fontsize=6.9,color='black' if value=='tf_s' else 'white' if v < -15 else 'black')
 fig.colorbar(im,ax=ax,shrink=.75,pad=.01)
finish(fig,'08_malla_suministro')
# F09 Span density errors at 300K
fig,axs=plt.subplots(1,2,figsize=(7.25,2.75));d=H[(H.origen=='malla')&(H.T_K==300)]
for m in models:
 g=d[d.modelo==m].sort_values('P_bar');axs[0].plot(g.P_bar,g.error_rho_pct,marker='o',ms=3.4,c=COL[m],label=m,lw=1.6);axs[1].plot(g.P_bar,g.error_h_J_kg/1000,marker='o',ms=3.4,c=COL[m],lw=1.6)
axs[0].axhline(0,color='0.4',lw=.7);axs[1].axhline(0,color='0.4',lw=.7)
axs[0].set(xlabel='Presión [bar]',ylabel='Error de densidad [%]',title='(a) Error frente a HEOS');axs[1].set(xlabel='Presión [bar]',ylabel='Error en entalpía residual [kJ/kg]',title='(b) Error de propiedad residual')
for ax in axs:styling(ax)
fig.legend(*axs[0].get_legend_handles_labels(),loc='lower center',bbox_to_anchor=(.5,-.085),ncol=3,frameon=False)
finish(fig,'09_contraste_heos')
# F10 RK4 convergence
fig,axs=plt.subplots(1,2,figsize=(7.25,2.75));df=RK
axs[0].loglog(df.paso_fijo_s,df.error_Tgas_K,'o-',color='#215a95',lw=1.8,label='RK4, error en T')
ref=df.error_Tgas_K.iloc[-2]*(df.paso_fijo_s/df.paso_fijo_s.iloc[-2])**4
axs[0].loglog(df.paso_fijo_s,ref,'--',color='#ba6334',label='Pendiente teórica 4')
axs[1].plot(df.paso_fijo_s[1:],df.orden_estimado_T[1:],'o-',color='#0b8379',lw=1.8)
axs[1].axhline(4,c='0.4',ls='--',lw=.9)
axs[0].set(xlabel='Paso fijo [s] (escala log)',ylabel='Error absoluto en T [K]',title='(a) Convergencia de RK4');axs[1].set(xlabel='Paso fijo [s]',ylabel='Orden observado',title='(b) Aproximación a orden 4')
for ax in axs:styling(ax)
axs[0].legend(frameon=False,loc='best',fontsize=7.5)
finish(fig,'10_convergencia')
# F11 energy partition
fig,axs=plt.subplots(1,2,figsize=(7.25,2.75));d=B[B.version=='V3'];xx=np.arange(3)
axs[0].bar(xx-.16,d.Qgw_kJ,width=.32,color='#0b8379',label='Gas → pared');axs[0].bar(xx+.16,d.Qwa_kJ,width=.32,color='#ba6334',label='Pared → ambiente')
axs[1].bar(xx,d.Qgw_kJ-d.Qwa_kJ,color='#215a95',width=.48)
axs[0].set(xlabel='Ecuación de estado',ylabel='Energía [kJ]',title='(a) Transferencias integradas');axs[1].set(xlabel='Ecuación de estado',ylabel='Energía acumulada en pared [kJ]',title='(b) Almacenamiento neto')
for ax in axs:ax.set_xticks(xx,labels=['Ideal','VdW','PR']);styling(ax)
axs[0].legend(frameon=False,fontsize=7.4)
finish(fig,'11_energia_pared')
# F12 iso limits
fig,axs=plt.subplots(1,2,figsize=(7.25,2.75));xx=np.arange(3);width=.34
axs[0].bar(xx-width/2,ISO.tf_adiabatico_s,width,label='Adiabático',color='#ba6334');axs[0].bar(xx+width/2,ISO.tf_isotermo_s,width,label='Isotermo',color='#215a95')
axs[1].bar(xx-width/2,ISO.mf_adiabatico_kg,width,color='#ba6334');axs[1].bar(xx+width/2,ISO.mf_isotermo_kg,width,color='#215a95')
for ax in axs:ax.set_xticks(xx,labels=['Ideal','VdW','PR']);styling(ax)
axs[0].set(ylabel='Tiempo [s]',title='(a) Duración');axs[1].set(ylabel='Masa final [kg]',title='(b) Masa retenida')
axs[0].legend(frameon=False,fontsize=7.4)
finish(fig,'12_isotermo_adiabatico')

print('Figuras generadas:',len(list(F.glob('*.pdf'))))
