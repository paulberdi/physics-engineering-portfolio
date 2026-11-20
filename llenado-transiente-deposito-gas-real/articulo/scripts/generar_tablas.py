from pathlib import Path
import pandas as pd
D=Path(__file__).resolve().parents[1]/'datos'
O=Path(__file__).resolve().parent/'tablas'; O.mkdir(exist_ok=True)
read=lambda n:pd.read_csv(D/n,encoding='utf-8-sig')
B=read('casos_base.csv'); C=read('comparacion_controlada_V1_V2.csv'); M=read('sensibilidad_masa_pared.csv'); S=read('sensibilidad_pared_3EOS.csv'); I=read('V2_limite_isotermo.csv'); H=read('propiedades_contra_HEOS.csv'); R=read('orden_RK4_paso_fijo_V3.csv');
rename={'Ideal':'Ideal','Van der Waals':'VdW','Peng–Robinson':'PR'}
def write(name,header,fmt,records):
 lines=[r'\begin{tabular}{'+fmt+'}',r'\toprule',header+r' \\',r'\midrule']
 for line in records: lines.append(' & '.join(str(x) for x in line)+r' \\')
 lines+=[r'\bottomrule',r'\end{tabular}']
 (O/name).write_text('\n'.join(lines)+'\n',encoding='utf8')
write('base.tex',r'Etapa & EOS & $P_f$ [bar] & $t_f$ [s] & $T_f$ [K] & $m_f$ [kg] & $Z_f$', 'llrrrrr',[[x.version,rename[x.modelo],f'{x.objetivo_bar:.0f}',f'{x.tf_s:.2f}',f'{x.Tf_K:.2f}',f'{x.mf_kg:.3f}',f'{x.Zf:.3f}'] for x in B.itertuples()])
write('control.tex',r'EOS & $P_f$ [bar] & V1 [s] & V2 [s] & $\Delta t$ [s] & $\Delta t/t_{V1}$ [\%]','lrrrrr',[[rename[x.eos],f'{x.Pobj_bar:.0f}',f'{x.t_V1_s:.2f}',f'{x.t_V2_s:.2f}',f'{x.diferencia_s:.2f}',f'{x.diferencia_pct:.2f}'] for x in C.itertuples()])
write('pared.tex',r'EOS & $m_w^*$ [kg] & $C_w$ [kJ/K] & $t_f$ [s] & $T_f$ [K] & $T_{w,f}$ [K]','lrrrrr',[[rename[x.eos],f'{x.masa_efectiva_kg:.0f}',f'{x.C_pared_J_K/1000:.1f}',f'{x.t_final_s:.2f}',f'{x.temperatura_final_K:.2f}',f'{x.temperatura_pared_final_K:.2f}'] for x in M.itertuples()])
write('transferencias.tex',r'EOS & $Q_{gw}$ [kJ] & $Q_{wa}$ [kJ] & $\Delta U_w$ [kJ] & $T_{w,f}$ [K]','lrrrr',[[rename[x.modelo],f'{x.Qgw_kJ:.2f}',f'{x.Qwa_kJ:.2f}',f'{x.Qgw_kJ-x.Qwa_kJ:.2f}',f'{x.Twf_K:.2f}'] for x in B[B.version=='V3'].itertuples()])
write('isotermo.tex',r'EOS & $t_{ad}$ [s] & $t_{iso}$ [s] & $m_{ad}$ [kg] & $m_{iso}$ [kg] & $Q_{ext}$ [kJ]','lrrrrr',[[rename[x.modelo],f'{x.tf_adiabatico_s:.2f}',f'{x.tf_isotermo_s:.2f}',f'{x.mf_adiabatico_kg:.3f}',f'{x.mf_isotermo_kg:.3f}',f'{x.calor_extraido_kJ:.1f}'] for x in I.itertuples()])
write('ua.tex',r'$UA_{gw}$ [W/K] & EOS & $t_f$ [s] & $T_{g,max}$ [K] & $Q_{gw}$ [kJ]','rlrrr',[[f'{x.valor:.0f}',rename[x.modelo],f'{x.tf_s:.2f}',f'{x.Tmax_K:.2f}',f'{x.Qgw_kJ:.1f}'] for x in S[S.parametro=='UA_gas_pared'].sort_values(['valor','modelo']).itertuples()])
write('heos.tex',r'EOS & $\rho$ [kg/m$^3$] & $e_\rho$ [\%] & $e_{h^r}$ [kJ/kg]','lrrr',[[rename[x.modelo],f'{x.rho_kg_m3:.2f}',f'{x.error_rho_pct:+.2f}',f'{x.error_h_J_kg/1000:+.2f}'] for x in H[(H.P_bar==200)&(H.T_K==300)].itertuples()])
write('rk.tex',r'$h$ [s] & Error $T$ [K] & Orden observado','rrr',[[f'{x.paso_fijo_s:.4f}',f'{x.error_Tgas_K:.3e}',f'{x.orden_estimado_T:.3f}' if pd.notna(x.orden_estimado_T) else '---'] for x in R.itertuples()])
print('Tablas',len(list(O.glob('*.tex'))))
