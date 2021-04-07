*****************************************************************
*
*                           READ ME
*
*         Copyright 2019, Tesfaye N. Asfaw, Leonid Tyan,
*           Alexey V. Glukhov, Vladimir E. Bondarenko
* Program m9da.f for mouse right atrial action potential model
*              Room temperature +25 C (298 K) control
*
*  Mathematical model consists of one file written in FORTRAN 90
*
*  Input data: initial conditions for model variables and model
*              parameters
*  Output data: action potential, ionic currents, model variables
*               as functions of time, simulated characterisctics
*
*  Program code consists of three modules: 
*  main module, subroutine 'RK4', subroutine 'fun'.
*  Subroutine 'RK4' is called from the main module to employ
*  4th order fixed step Runge-Kutta method.
*  Subroutine 'fun' is called from 'RK4' and describes right-hand
*  parts of ordinary differential equations and model parameters. 
*
*  To compile a code use command:
*
*         ifort -o m9da -O3 -i4 -r8 m9da.f
*
*  Code is executed in the same directory as:
*
*                   ./m9da
*
*****************************************************************
      program mouse
      implicit none

      integer NEQ,NT,i,ii,iii,j,NP,IP
      real time,Vrest,fun,TINI,TPER,TW,AMPL,STON(5000),STOFF(5000)
      real dt,durat,stimon1,stimoff1,stimst1,stimon2,stimoff2,stimst2
      real sum_i,Istim,Icas,Icab,Inaca,Ipca,Ina,Inab,Iclca,Icasc,Icase
      real Inak,Ikto,Ik1,Iks,Ikur,Ikr,Ikur1,Ikur2,Ikur3,Ikur4
      real Jrelt,Jtrt,Jxfert,Jleak,Jup,Jtrpn
      real Acap,Vmyo,Vjsr,Vnsr,Vss,Cm,F,Bi,Bss,Bjsr
      real frel,fcal,fup,ftrpn,fnaca,fxfer,fleak,fcab,fpca,frelm,fcalm
      real srel,scal,sup,strpn,snaca,sxfer,sleak,scab,spca,gain(30)
      real fnaf,fnab,fnanaca,fnanak,fkv,fknak,fkstim,fcat
      real snaf,snab,snanaca,snanak,skv,sknak,skstim,scat
      real sikto,sik1,siksikr,sikur2,sikur3
      real Y(300),YDOT(300),Ts(30000),Inad(30,30000),Y2ini
      real Vm(30),Inasm1(30),Inasm2(30),Inasn1(30),Inasn2(30)
      real Inasnm1,Inasnm2,Inast1,Inast2,Caim(30),Cain(30),Cainm,Caism
      real G1(30),G1n(30),G1nm
      real RCcavf,RCecavf,RCcytf,cAMPcell,Ccell,ACactcell,PDEactcell
      real Gscavabg,Gsecavabg,Gscytabg,RB1cavtot,RB1ecavtot
      real RB1cyttot,L,dL
      real Gicavabg,Giecavabg,RB2cavtot,RB2ecavtot
      real RB2pPKA,RB2pbARK,RB2ptot
      real PDE2cavtot,PDE3cavtot,PDE4cavtot,PDE2ecavtot,
     *  PDE3ecavtot,PDE4ecavtot,PDE2cyttot,PDE3cyttot,PDE4cyttot,
     *  PDE8cavtot,PDE8ecavtot,PDE8cyttot
      real fluxecavcyt,fluxcytcav,fluxcytecav
      real fluxcavecav,fluxcavcyt,fluxecavcav
      real PDEact2,PDEact3,PDEact4,PDEact8,PKAactcell
      real PKAcav,PKAecav,PKAcyt,PKIcavf,PKIecavf,PKIcytf
      real RB1pPKA,RB1pbARK,RB1ptot,YICaLcc,YICaLcp,YRyRc,YRyRp
      real YINac,YINap,pp1cytf,YICaLec,YICaLep,YICaLc,YICaLp,fICaLcav
      real PDE2mem,PDE3mem,PDE4mem,PDE8mem,PDEmem
      real RB2GicavPKA,LRB2GicavPKA,RB2cavPKAf,Gicavf,
     * acavB2i,bcavB2i,ccavB2i
      real RB2GiecavPKA,LRB2GiecavPKA,RB2ecavPKAf,Giecavf,
     * aecavB2i,becavB2i,cecavB2i
      real Icat,Icatc,Icats,ICaK
      
      EXTERNAL FUN
      
      common/a1/Jrelt,Jtrt,Jxfert,Jleak,Jup,Jtrpn
      common/a2/Acap,Vmyo,Vjsr,Vnsr,Vss,Cm,F,Bi,Bss,Bjsr
      common/a10/Icas,Icab,Inaca,Ipca,Ina,Inab,Iclca,Icasc,Icase
      common/a11/Inak,Ikto,Ik1,Iks,Ikur,Ikr,Ikur1,Ikur2,Ikur3,Ikur4
      common/b1/RCcavf,RCecavf,RCcytf,Gscavabg,Gsecavabg,Gscytabg
      common/b2/cAMPcell,Ccell,ACactcell,PDEactcell
      common/b3/fluxcavecav,fluxcavcyt,fluxecavcav              
      common/b4/fluxecavcyt,fluxcytcav,fluxcytecav 
      common/b5/RB1cavtot,RB1ecavtot,RB1cyttot
      common/b6/PDE2cavtot,PDE3cavtot,PDE4cavtot,PDE2ecavtot,
     *  PDE3ecavtot,PDE4ecavtot,PDE2cyttot,PDE3cyttot,PDE4cyttot,
     *  PDE8cavtot,PDE8ecavtot,PDE8cyttot
      common/b7/PDEact2,PDEact3,PDEact4,PDEact8,PKAactcell
      common/b8/PKAcav,PKAecav,PKAcyt,PKIcavf,PKIecavf,PKIcytf
      common/b9/RB1pPKA,RB1pbARK,RB1ptot,L,RB2pPKA,RB2pbARK,RB2ptot
      common/b10/PDE2mem,PDE3mem,PDE4mem,PDE8mem,PDEmem
      common/b11/pp1cytf,fICaLcav
      common/b12/Gicavabg,Giecavabg,RB2cavtot,RB2ecavtot
      common/c1/RB2GicavPKA,LRB2GicavPKA,RB2cavPKAf,Gicavf,
     * acavB2i,bcavB2i,ccavB2i
      common/c2/RB2GiecavPKA,LRB2GiecavPKA,RB2ecavPKAf,Giecavf,
     * aecavB2i,becavB2i,cecavB2i
      common/d1/Icat,Icatc,Icats,ICaK

      NEQ = 185
      NT  = 13
     
      dL = sqrt(sqrt(10.0))
      L = 1.0E-04/dL
      do 233 j=1,1
      L = L*dL
      L = 0.0
                        
      dt       =     0.0001
      iii      =  1000
      ii       =     0
      durat    =300000.0
      stimon1  =    20.0
      stimoff1 =    20.5
      stimst1  =    75.0
      stimon2  =   520.0
      stimoff2 =  1020.0
      stimst2  =    20.0
      Vrest    =   -80.0
*
          NP   = 2000
          TINI = 20
          TPER = 1000
          TW   = 1.0
          AMPL = 80.0
          DO i=1,NP
          STON(i)=TINI+(i-1)*TPER
          STOFF(i)=TINI+(i-1)*TPER+TW
          END DO
*
* Initial Conditions
*
      frel  = 0.0
      fcal  = 0.0
      fup   = 0.0
      ftrpn = 0.0
      fnaca = 0.0
      fxfer = 0.0
      fleak = 0.0
      fcab  = 0.0
      fpca  = 0.0
      frelm = 0.0
      fcalm = 0.0
      fnaf  = 0.0
      fnab  = 0.0
      fnanaca = 0.0
      fnanak = 0.0
      fkv   = 0.0
      fknak = 0.0
      fkstim = 0.0
*
      srel  = 0.0
      scal  = 0.0
      sup   = 0.0
      strpn = 0.0
      snaca = 0.0
      sxfer = 0.0
      sleak = 0.0
      scab  = 0.0
      spca  = 0.0
      snaf  = 0.0
      snab  = 0.0
      snanaca = 0.0
      snanak = 0.0
      skv   = 0.0
      sknak = 0.0
      skstim = 0.0
      sikto =0.0
      sik1 = 0.0
      siksikr = 0.0
      sikur2 = 0.0
      sikur3 = 0.0
*
      time  =      0.0
      Y(1)  =  -0.8088310412E+02
      Y(2)  =   0.9279026900E-01
      Y(3)  =   0.1620528775E+04
      Y(4)  =   0.8053999501E+01
      Y(5)  =   0.1222160781E+03
      Y(6)  =   0.1620528775E+04
      Y(7)  =   0.9279026946E-01
      Y(8)  =   0.8851691228E-12
      Y(9)  =   0.9615600724E+00
      Y(10) =   0.3767463428E-02
      Y(11) =   0.5535450499E-05
      Y(12) =   0.3614717015E-08
      Y(13) =   0.7908511386E-12
      Y(14) =   0.7153680173E-08
      Y(15) =   0.6391662360E-08
      Y(16) =   0.9967529680E+00
      Y(17) =   0.7087660992E-04
      Y(18) =   0.6300248702E-05
      Y(19) =   0.2112480482E-10
      Y(20) =   0.4628585970E+00
      Y(21) =   0.1062737990E-01
      Y(22) =   0.9600858225E-04
      Y(23) =   0.1376322742E+05
      Y(24) =   0.1418267743E+06
      Y(25) =   0.4139033547E-02
      Y(26) =   0.9999623535E+00
      Y(27) =   0.3298367341E-03
      Y(28) =   0.5091689794E-03
      Y(29) =   0.9980927689E+00
      Y(30) =   0.9979010383E+00
      Y(31) =   0.1113258849E-02
      Y(32) =   0.7192382272E-03
      Y(33) =   0.2223437325E-03
      Y(34) =   0.4412087393E-04
      Y(35) =   0.5091689794E-03
      Y(36) =   0.9980927689E+00
      Y(37) =   0.1552405441E-06
      Y(38) =   0.6461071366E-04
      Y(39) =   0.3762774987E-05
      Y(40) =   0.8212756509E-08
      Y(41) =   0.7151890060E-02
      Y(42) =   0.3114893963E+00
      Y(43) =   0.5091689794E-03
      Y(44) =   0.3779115232E-12
      Y(45) =   0.1193078173E-02
      Y(46) =   0.3630488755E-61
      Y(47) =   0.3490326981E-02
      Y(48) =   0.5883268221E-02
      Y(49) =   0.3955088510E-03
      Y(50) =   0.4450890836E-01
      Y(51) =   0.3630488755E-61
      Y(52) =   0.2392020281E-01
      Y(53) =   0.6409630097E-02
      Y(54) =   0.2487944176E-02
      Y(55) =   0.1722672905E-02
      Y(56) =   0.3630488755E-61
      Y(57) =   0.3771728607E-03
      Y(58) =   0.7236469606E-03
      Y(59) =   0.3474740964E-03
      Y(60) =   0.1300635182E+06
      Y(61) =   0.5077713406E+05
      Y(62) =   0.2353396779E+05
      Y(63) =   0.2193421483E+04
      Y(64) =   0.1951294807E-01
      Y(65) =   0.9059006924E-02
      Y(66) =   0.1097267542E+05
      Y(67) =   0.9355062972E+05
      Y(68) =   0.2873599432E+05
      Y(69) =   0.6773138136E-90
      Y(70) =   0.1449457584E-01
      Y(71) =   0.1340498892E+05
      Y(72) =   0.2242397227E+00
      Y(73) =   0.3690967042E+05
      Y(74) =   0.1340729197E-02
      Y(75) =   0.4134187505E-02
      Y(76) =   0.4430125145E+04
      Y(77) =   0.7006267909E+04
      Y(78) =   0.1421393882E+05
      Y(79) =   0.7736146378E+01
      Y(80) =   0.3812957050E+00
      Y(81) =   0.5158797692E-01
      Y(82) =   0.8896934822E+00
      Y(83) =   0.7537918546E-01
      Y(84) =   0.8253133062E+00
      Y(85) =   0.6836347022E+01
      Y(86) =   0.6205890724E+00
      Y(87) =   0.1172793205E+00
      Y(88) =   0.1154255173E+01
      Y(89) =   0.1320878782E+00
      Y(90) =   0.1033166329E+01
      Y(91) =   0.9274466720E+01
      Y(92) =   0.1086997100E+00
      Y(93) =   0.1546752825E-01
      Y(94) =   0.2929517738E+00
      Y(95) =   0.6863855598E-01
      Y(96) =   0.2353122184E+00
      Y(97) =   0.3382412669E+00
      Y(98) =   0.4724516017E+00
      Y(99) =   0.4126582483E+00
      Y(100) =  0.2320312876E-01
      Y(101) =  0.6026039427E+00
      Y(102) =  0.6026039427E+00
      Y(103) =  0.2596146734E+00
      Y(104) =  0.3675832948E+00
      Y(105) =  0.3813640007E+00
      Y(106) =  0.8851705616E-12
      Y(107) =  0.2605069110E-10
      Y(108) =  0.3415627435E-01
      Y(109) =  0.5076951530E-03
      Y(110) =  0.2829869165E-05
      Y(111) =  0.7010470517E-08
      Y(112) =  0.2327573109E-10
      Y(113) =  0.5549734835E-07
      Y(114) =  0.4958553899E-07
      Y(115) =  0.6512671925E-11
      Y(116) =  0.3169194559E-02
      Y(117) =  0.6261242517E-06
      Y(118) =  0.3339118148E-07
      Y(119) =  0.1866017795E-12
      Y(120) =  0.1213360465E+00
      Y(121) =  0.2785949016E-02
      Y(122) =  0.2516918244E-04
      Y(123) =  0.4069800161E-07
      Y(124) =  0.1694073185E-04
      Y(125) =  0.1028703688E-05
      Y(126) =  0.6604024237E-08
      Y(127) =  0.1875148561E-02
      Y(128) =  0.8166786118E-01
      Y(129) =  0.3801202195E+00
      Y(130) =  0.5091689794E-03
      Y(131) =  0.9980927689E+00
      Y(132) =  0.9214774521E+00
      Y(133) =  0.9908502903E+00
      Y(134) =  0.2087704319E+00
      Y(135) =  0.8637307739E-03
      Y(136) =  0.9999881755E+00
      Y(137) =  0.8187695769E-12
      Y(138) =  0.8894320487E+00
      Y(139) =  0.3484859972E-02
      Y(140) =  0.5120227116E-05
      Y(141) =  0.3343570256E-08
      Y(142) =  0.8187721513E-12
      Y(143) =  0.7315039113E-12
      Y(144) =  0.6616944790E-08
      Y(145) =  0.5912117549E-08
      Y(146) =  0.8046441685E-10
      Y(147) =  0.1055007176E+00
      Y(148) =  0.1568150421E-02
      Y(149) =  0.8740796912E-05
      Y(150) =  0.2165368505E-07
      Y(151) =  0.2011610269E-10
      Y(152) =  0.7189314575E-10
      Y(153) =  0.1714179939E-06
      Y(154) =  0.1531578593E-06
      Y(155) =  0.6077748720E-01
      Y(156) =  0.5796341475E-38
      Y(157) =  0.3591463241E-02
      Y(158) =  0.4069691772E-03
      Y(159) =  0.4580528431E-03
      Y(160) =  0.5796341475E-39
      Y(161) =  0.2305390265E-05
      Y(162) =  0.2397840156E-06
      Y(163) =  0.8422801321E+00
      Y(164) =  0.5090725598E-01
      Y(165) =  0.1153809439E-02
      Y(166) =  0.1162267165E-04
      Y(167) =  0.4390450831E-07
      Y(168) =  0.6508479153E-06
      Y(169) =  0.6016286658E-01
      Y(170) =  0.3636232570E-01
      Y(171) =  0.8241495995E-02
      Y(172) =  0.8301908322E-03
      Y(173) =  0.3136036308E-05
      Y(174) =  0.4648913681E-04
      Y(175) =  0.2965524523E-01
      Y(176) =  0.2144114378E-01
      Y(177) =  0.8608445698E+00
      Y(178) =  0.1198169988E+00
      Y(179) =  0.1334142186E-01
      Y(180) =  0.2971089895E-03
      Y(181) =  0.2134627497E-02
      Y(182) =  0.3565307875E-02
      Y(183) =  0.1184324373E-01
      Y(184) =  0.5091689794E-03
      Y(185) =  0.9980927689E+00
      Y2ini =   Y(2)
*
*  Read initial conditions from file
*
*      OPEN(UNIT=101,FILE='fort.101',STATUS='OLD')
*
*      do i=1,NEQ
*      read(101,*) Y(i)
*      write(103,45) Y(i)
*      end do
*
*      CLOSE(UNIT=101)
*
      do 65 i=1,NEQ
      YDOT(i) = 0.0
  65  continue
*
  77  continue
      if (time.ge.durat) go to 88
      if(mod(ii,iii).ne.0) go to 6
*      if(time.le.95000.0) go to 6
*      write(6,44) time,Y(1),Y(2),Y(3),Y(4),Y(5),Y(6)
      YICaLcc = Y(8)+Y(9)+Y(10)+Y(11)+Y(12)+Y(13)+Y(14)+Y(15)+Y(106)
      YICaLcp = Y(107)+Y(108)+Y(109)+Y(110)+Y(111)+Y(112)+Y(113)
     *       + Y(114)+Y(115)
      YICaLec = Y(137)+Y(138)+Y(139)+Y(140)+Y(141)+Y(142)+Y(143)
     *       + Y(144)+Y(145)
      YICaLep = Y(146)+Y(147)+Y(148)+Y(149)+Y(150)+Y(151)+Y(152)
     *       + Y(153)+Y(154)
      YICaLc = fICaLcav*YICaLcc+(1-fICaLcav)*YICaLec
      YICaLp = fICaLcav*YICaLcp+(1-fICaLcav)*YICaLep
      YRyRc  = Y(16)+Y(17)+Y(18)+Y(19)
      YRyRp  = Y(116)+Y(117)+Y(118)+Y(119)
      YINac  = Y(20)+Y(21)+Y(22)+Y(37)+Y(38)+Y(39)+Y(40)+Y(41)+Y(42)
      YINap  = Y(120)+Y(121)+Y(122)+Y(123)+Y(124)+Y(125)+Y(126)
     *       + Y(127)+Y(128)
      write(11,44) time,Y(1),Icas,Icab,Inaca,Ipca,Ina,Inab,Icasc,Icase
      write(12,44) time,Y(1),Inak,Ikto,Ik1,Iks,Ikur,Ikr,Iclca,Istim,
     *  Icat,Icatc,Icats,ICaK
      write(13,44) time,Y(1),Y(2),Y(3),Y(4),Y(5),Y(6)
      write(14,44) time,Y(7),Y(8),Y(9),Y(10),Y(11),Y(12)
      write(15,44) time,Y(13),Y(14),Y(15),Y(16),Y(17),Y(18)
      write(16,44) time,Y(19),Y(20),Y(21),Y(22),Y(23),Y(24)
      write(17,44) time,Y(25),Y(26),Y(27),Y(28),Y(29),Y(30)
      write(18,44) time,Y(31),Y(32),Y(33),Y(34),Y(35),Y(36)
      write(19,44) time,Y(37),Y(38),Y(39),Y(40),Y(41),Y(42)
      write(10,44) time,Y(43),Y(44),Y(18)+Y(19)
      write(33,44) time,Y(45),Y(46),Y(47),Y(48),Y(49),Y(50)
      write(34,44) time,Y(51),Y(52),Y(53),Y(54),Y(55),Y(56)
      write(35,44) time,Y(57),Y(58),Y(59),Y(60),Y(61),Y(62)
      write(36,44) time,Y(63),Y(64),Y(65),Y(66),Y(67),Y(68)
      write(37,44) time,Y(69),Y(70),Y(71),Y(72),Y(73),Y(74)
      write(38,44) time,Y(75),Y(76),Y(77),Y(78),Y(79),Y(80)
      write(39,44) time,Y(81),Y(82),Y(83),Y(84),Y(85),Y(86)
      write(40,44) time,Y(87),Y(88),Y(89),Y(90),Y(91),Y(92)
      write(41,44) time,Y(93),Y(94),Y(95),Y(96),Y(97),Y(98)
      write(42,44) time,Y(99),Y(100),Y(101),Y(102),Y(103),Y(104)
      write(43,44) time,Y(105),Y(106),Y(107),Y(108),Y(109),Y(110)
      write(44,44) time,Y(111),Y(112),Y(113),Y(114),Y(115),Y(116)
      write(45,44) time,Y(117),Y(118),Y(119),Y(120),Y(121),Y(122)
      write(46,44) time,Y(123),Y(124),Y(125),Y(126),Y(127),Y(128)
      write(47,44) time,Y(129),Y(130),Y(131),Y(132),Y(133),Y(134)
      write(48,44) time,Y(135),Y(136),Y(137),Y(138),Y(139),Y(140)
      write(49,44) time,Y(141),Y(142),Y(143),Y(144),Y(145),Y(146)
      write(62,44) time,Y(147),Y(148),Y(149),Y(150),Y(151),Y(152),
     * Y(153),Y(154)
      write(64,44) time,Y(155),Y(156),Y(157),Y(158),Y(159),Y(160),
     * Y(161),Y(162)
      write(67,44) time,Y(163),Y(164),Y(165),Y(166),Y(167),Y(168),
     * Y(169),Y(170),Y(171),Y(172),Y(173),Y(174),Y(175),Y(176)
      write(68,44) time,Y(177),Y(178),Y(179),Y(180),Y(181),Y(182),
     * Y(183),Y(184),Y(185)
      write(50,44) time,RCcavf,RCecavf,RCcytf,cAMPcell,Ccell,ACactcell,
     *   PDEactcell,PKAactcell,pp1cytf
      write(20,44) time,Y(1),Ikto+Iks+Ikur+Ikr,Ikur1,Ikur2,Ikur3,Ikur4
      write(30,44) time,frel,fcal,fup,ftrpn,fnaca,fxfer,fleak,fcab,fpca,
     *  fnaf,fnab,fnanaca,fnanak,fkv,fknak,fkstim,fcat
      write(31,44) time,Bi,Bss,Bjsr
      write(32,44) time,srel,scal,sup,strpn,snaca,sxfer,sleak,scab,spca,
     *  snaf,snab,snanaca,snanak,skv,sknak,skstim,
     *  sikto,sik1,siksikr,sikur2,sikur3,scat
      write(51,44) time,fluxcavecav,fluxcavcyt,fluxecavcav 
     *             ,fluxecavcyt,fluxcytcav,fluxcytecav 
      write(52,44) time,YDOT(60),YDOT(61),YDOT(62),YDOT(63)
      write(53,44) time,YDOT(66),YDOT(67),YDOT(68),YDOT(71),YDOT(72),
     *   YDOT(73),YDOT(76),YDOT(77),YDOT(78),YDOT(101),YDOT(102),
     *   YDOT(103)
      write(54,44) time,YDOT(79),YDOT(85),YDOT(91)
      write(55,44) time,Gscavabg,Gsecavabg,Gscytabg,RB1cavtot,
     *   RB1ecavtot,RB1cyttot,RB1pPKA,RB1pbARK,RB1ptot,Gicavabg,
     *   Giecavabg,RB2cavtot,RB2ecavtot,RB2pPKA,RB2pbARK,RB2ptot
      write(56,44) time,PDE2cavtot,PDE3cavtot,PDE4cavtot,PDE2ecavtot,
     *  PDE3ecavtot,PDE4ecavtot,PDE2cyttot,PDE3cyttot,PDE4cyttot,
     *  PDE8cavtot,PDE8ecavtot,PDE8cyttot
      write(57,44) time,PDEact2,PDEact3,PDEact4,PDEact8
      write(58,44) time,PKAcav,PKAecav,PKAcyt,PKIcavf,PKIecavf,PKIcytf
      write(59,44) time,PDE2mem,PDE3mem,PDE4mem,PDE8mem,PDEmem
      write(60,44) time,YICaLc,YICaLp,YRyRc,YRyRp,YINac,YINap
      write(61,44) time,Y(133),(1-Y(133)),Y(134),(1-Y(134))
      write(63,44) time,YICaLcc,YICaLcp,YICaLec,YICaLep
      write(65,44) time,RB2GicavPKA,LRB2GicavPKA,RB2cavPKAf,Gicavf,
     * Gicavabg,acavB2i,bcavB2i,ccavB2i
      write(66,44) time,RB2GiecavPKA,LRB2GiecavPKA,RB2ecavPKAf,
     * Giecavf,Giecavabg,aecavB2i,becavB2i,cecavB2i
*
   6  continue
      ii = ii + 1
*
        DO I=1,NP
        IF((time.GE.STON(I)).AND.(time.LE.STOFF(I))) THEN
              Istim = AMPL
              IP = I
              GO TO 78
        ELSE
              Istim = 0.0
        END IF
        END DO
  78    CONTINUE
        DO I=1,NP
        IF((time.GE.STON(I)).AND.(time.LE.STON(I)+15.0)) THEN
              dt  = 0.000002
              iii = 250000
              GO TO 79
        ELSE
              dt  = 0.0001
              iii = 5000
        END IF
        END DO
   79   CONTINUE
*
       Y(9)  = 1.0-(Y(8)+Y(10)+Y(11)+Y(12)+Y(13)+Y(14)+Y(15)+Y(106)
     *       + Y(107)+Y(108)+Y(109)+Y(110)+Y(111)+Y(112)+Y(113)
     *       + Y(114)+Y(115))
       Y(138)= 1.0-(Y(137)+Y(139)+Y(140)+Y(141)+Y(142)+Y(143)+Y(144)
     *       + Y(145)+Y(146)+Y(147)+Y(148)+Y(149)+Y(150)+Y(151)+Y(152)
     *       + Y(153)+Y(154))
*       Y(16) = 1.0-(Y(17)+Y(18)+Y(19)+Y(116)+Y(117)+Y(118)+Y(119))
       Y(20) = 1.0-(Y(21)+Y(22)+Y(37)+Y(38)+Y(39)+Y(40)+Y(41)+Y(42)
     *       + Y(120)+Y(121)+Y(122)+Y(123)+Y(124)+Y(125)+Y(126)
     *       + Y(127)+Y(128))
       Y(30) = 1.0-(Y(31)+Y(32)+Y(33)+Y(34))
*
       call RK4(FUN,NEQ,Y,YDOT,TIME,dt,Istim)
*
      frel  = Jrelt*Vjsr/Vmyo
      fcal  = -1.0*Icas*Acap*Cm/(2.0*F*Vmyo)
      fcat  = -1.0*Icat*Acap*Cm/(2.0*F*Vmyo)
      fup   = -1.0*Jup
      ftrpn = -1.0*Jtrpn
      fnaca = 2.0*Inaca*Acap*Cm/(2.0*F*Vmyo)
      fxfer = Jxfert
      fleak = Jleak
      fcab  = -1.0*Icab*Acap*Cm/(2.0*F*Vmyo)
      fpca  = -1.0*Ipca*Acap*Cm/(2.0*F*Vmyo)
*
      fnaf  = -1.0*Ina*Acap*Cm/(F*Vmyo)
      fnab  = -1.0*Inab*Acap*Cm/(F*Vmyo)
      fnanaca = -1.0*3.0*Inaca*Acap*Cm/(F*Vmyo)
      fnanak = -1.0*3.0*Inak*Acap*Cm/(F*Vmyo)
*
      fkv = -1.0*(Ikto+Ik1+Iks+Ikur+Ikr+ICaK)*Acap*Cm/(F*Vmyo)
      fknak = -1.0*(-2.0)*Inak*Acap*Cm/(F*Vmyo)
      fkstim = -1.0*(-1.0)*Istim*Acap*Cm/(F*Vmyo)
*
        if((time.ge.299001).and.(time.le.300001)) then
      srel  = srel+frel*dt
      scal  = scal+fcal*dt
      sup   = sup+fup*dt
      strpn = strpn+ftrpn*dt
      snaca = snaca+fnaca*dt
      sxfer = sxfer+fxfer*dt
      sleak = sleak+fleak*dt
      scab  = scab+fcab*dt
      spca  = spca+fpca*dt
      scat  = scat+fcat*dt
*
      snaf  = snaf+fnaf*dt
      snab  = snab+fnab*dt
      snanaca = snanaca+fnanaca*dt
      snanak = snanak+fnanak*dt
*
      skv   = skv+fkv*dt
      sknak = sknak+fknak*dt
      skstim = skstim+fkstim*dt
      sikto = sikto+dt*(-1.0)*Ikto*Acap*Cm/(F*Vmyo)
      sik1 = sik1+dt*(-1.0)*Ik1*Acap*Cm/(F*Vmyo)
      siksikr = siksikr+dt*(-1.0)*(Iks+Ikr)*Acap*Cm/(F*Vmyo)
      sikur2 = sikur2+dt*(-1.0)*Ikur2*Acap*Cm/(F*Vmyo)
      sikur3 = sikur3+dt*(-1.0)*Ikur3*Acap*Cm/(F*Vmyo)
*
        end if
*
      go to 77
*
  88  continue
*
      write(150,44)L,YDOT(60),YDOT(61),YDOT(62),YDOT(63),
     *  RB1ptot,PDEact2,PDEact3,PDEact4,PDEact8,cAMPcell,Ccell,
     *  ACactcell,PDEactcell,PKAactcell,Y(104),Y(3),Y(105),
     *  Y(129),Inak,Ikur2,Ik1,Ikto,YICaLc,YICaLp,YRyRc,YRyRp,
     *  Y(89),Y(16),Y(116)
  233 continue
*
*   Write new initial conditions
*
      do j=1,NEQ
      write(102,45) Y(j)
      end do
*
  44  format(f14.5,40e15.6)
  45  format(30e20.10)
* END of Program!
      stop
      end
*****************************************************************
*
      subroutine fun(NEQ,time,Y,YDOT,Istim)
*
*****************************************************************
      implicit none
      integer NEQ,i
      real V,Va,Vi,time
      real temp,temp1,temp2,temp3,temp4,temp5,temp6,temp7
      real temp8,temp9,temp10,temp11,temp12,temp13
      real tempa1,tempa2,tempa3,tempa41,tempa42,tempa4
      real tempa5,tempa6,tempa7,tempa8,tempa9,tempa10,tempa11
      real Bi,Bss,Bjsr,sigma
      real alp11,bet11,alp12,bet12,alp13,bet13
      real alp2,bet2,alp3,bet3,alp4,bet4,alp5,bet5
      real alp25,bet25,alp26,bet26,alp27,bet27
      real alp25p,bet25p,alp26p,bet26p,temp7p,temp8p,temp9p,temp10p
      real ass,iss1,taua1,ala0,bea0,ala1,bea1,ali,bei
      real iss2,taua2,taua3
      real Acap,Vcell,Vmyo,Vjsr,Vnsr,Vss,Vcav,Vecav,cKo,cNao,cCao
      real v1,v2,v3,Kmup,Ttr,Txfer,Kap,Kam,Kbp,Kbm,Kcp,Kcm
      real Kpcmax,Kpchalf,Kpcb,Kpcf,Vcyt
      real SFC4I1,SFC4I2,SFC4I3,SFOI1,SFOI2,SFI1I3,SFI2I3,SFICA
      real cLTRPNtot,cHTRPNtot,khtrpnp,khtrpnm,kltrpnp,kltrpnm
      real cCMDNtot,cCSQNtot,Kmcmdn,Kmcsqn
      real Cm,F,T,R,factor,ifactor,Gna,Gkp,Pnak,knaca,Kmna,Kmca,ksat
      real nu,Inakmax,Kmnai,Kmko,Ipcamax,Kmpca,Gcab,Gnab
      real Gkto,Gks,Gkur,kf,kb,Gkr,Gk1,alpha,beta,gammac,gammae
      real Gclca,Poclcamax,Kmclca,Ecl,Ecan,ENa,EK,EKr
      real sum_i,Istim,Icas,Icab,Inaca,Ipca,Ina,Inab,Iclca,Icasc,Icase
      real Inak,Ikto,Ik1,Iks,Ikur,Ikr,Ikur1,Ikur2,Ikur3,Ikur4
      real Jrelt,Jtrt,Jxfert,Jleak,Jup,Jtrpn,t1,t2,Icasmax
      real Gkur1,Gkur2,Gkur3,Gkur4,taui1,taui2,taui3,taui4
      real L,RB1tot,fcavB1,fecavB1,fcytB1,Gstot,fcavGs,fecavGs,fcytGs
      real KB1L,KB1H,KB1C,kPKAp,kPKAm,kBARKp,kBARKm,kact1Gs,kact2Gs
      real khydGs,kreasGs,Kmatp,ATP,ACtot,fac56ac47,fcavac56,factGsGi
      real RB2tot,fcavB2,fecavB2,Gitot,fcavGi,fecavGi
      real KB2L,KB2H,KB2C,KB2F,KB2A,kact1Gi,kact2Gi
      real Kac56mgsgi,Kac56mgi,hac56gsgi,Vac56gsgi
      real fecavac47,Kac56mgsa,hac56gsa,Vac56gbg,Kac56mgsbg,hac56gsbg
      real AC56bas,AF56,Kac47mgsa,hac47gsa,Vac47gbg,Kac47mgsbg
      real hac47gsbg,AC47bas,AF47,IBMX,hibmxpde2,Kibmxpde2
      real hibmxpde3,Kibmxpde3,hibmxpde4,Kibmxpde4,kfpdep,kbpdep
      real deltakpde34,kpde2,Kmpde2,kpde3,Kmpde3,kpde4,Kmpde4
      real fpdepart,rpartpde23,rpartpde34,PDE2tot,PDE3tot,PDE4tot
      real fcavpde2,fecavpde2,fcytpde2,fcavpde3,fecavpde3,fcytpde3
      real fcavpde4,fecavpde4,fcytpde4,PKAtot,fcavpka,fecavpka
      real fcytpka,PKItot,fcavpki,fecavpki,fcytpki,kpkaif1,kpkai1
      real kpkaif2,kpkai2,kpkaif3,kpkai3,kpkaiif1,kpkaii1,kpkaiif2
      real kpkaii2,kpkaiif3,kpkaii3,kpkif,kpki,pp1cyttot,inhib1cyttot
      real Kinhib1,kpkainhib1,Kmpkainhib1,kpp2ainhib1pp2acyt
      real Kmpp2ainhib1,RB1cavtot,Gscavabg,RB1cavnptot,acavB1,bcavB1
      real ccavB1,RB1cavnpf,Gscavf,LRB1cavnp,RB1Gscavnp,LRB1Gscavnp
      real RB2cavtot,Gicavabg,RB2cavnptot,acavB2i,bcavB2i,ccavB2i
      real RB2cavPKAf,Gicavf,LRB2cavPKA,RB2GicavPKA
      real LRB2GicavPKA,acavB2s,bcavB2s,ccavB2s,dcavB2s
      real pcav,qcav,rcav,Acav,Bcav,Dcav,Mcav,Ncav,z1cav,z2cav,z3cav
      real y1cav,y2cav,y3cav,ficav,pi
      real RB2cavnpf,LRB2cavnp,RB2Gscavnp,LRB2Gscavnp
      real RB2ecavtot,khydGi,kreasGi,RB2pPKA,RB2pbARK,RB2ptot
      real AC56cav,kcavAC56,PDE2cavtot,PDE3cavtot,PDE4cavtot
      real Kibmxmpde2,Kibmxmpde3,Kibmxmpde4,PKAcav,RCcavf,PKIcavf
      real kpkaiib1,kpkaiib2,kpkaiib3,kpkib,RB1ecavtot,Gsecavabg
      real RB1ecavnptot,aecavB1,becavB1,cecavB1,RB1ecavnpf,Gsecavf
      real LRB1ecavnp,RB1Gsecavnp,LRB1Gsecavnp,AC47ecav,kecavAC47
      real RB2ecavnptot,Giecavabg,aecavB2i,becavB2i,cecavB2i
      real RB2ecavPKAf,Giecavf,RB2ecavPKAtot,LRB2ecavPKA,RB2GiecavPKA
      real LRB2GiecavPKA,aecavB2s,becavB2s,cecavB2s,decavB2s
      real pecav,qecav,recav,Aecav,Becav,Decav,Mecav,Necav
      real y1ecav,y2ecav,y3ecav,fiecav,z1ecav,z2ecav,z3ecav
      real RB2ecavnpf,LRB2ecavnp,RB2Gsecavnp,LRB2Gsecavnp
      real PDE2ecavtot,PDE3ecavtot,PDE4ecavtot,PKAecav,RCecavf
      real PKIecavf,RB1cyttot,Gscytabg,RB1cytnptot,acytB1
      real bcytB1,ccytB1,RB1cytnpf,Gscytf,LRB1cytnp,RB1Gscytnp
      real LRB1Gscytnp,AC56cyt,kcytAC56,AC47cyt,kcytAC47,PDE2cyttot
      real PDE3cyttot,PDE4cyttot,PKAcyt,RCcytf,PKIcytf,kpkaib1,kpkaib2
      real kpkaib3,inhib1cytf,ainh1,binh1,Kinh1,cinh1,inhib1cytp
      real pp1cytf,jcavecav,jcavcyt,jecavcyt,cAMPcell,Ccell,ACactcell
      real fluxcavecav,fluxcavcyt,fluxecavcav,fluxecavcyt,fluxcytcav
      real fluxcytecav,PDEactcell,PDEact2,PDEact3,PDEact4,PKAactcell
      real RB1pPKA,RB1pbARK,RB1ptot,PDEact8
      real kpde8,Kmpde8,PDE8tot,fcavpde8,fecavpde8,fcytpde8
      real PDE8cavtot,PDE8ecavtot,PDE8cyttot
      real PDE2mem,PDE3mem,PDE4mem,PDE8mem,PDEmem
      real kPLBPKA,KPLBmPKA,kPLBPP1,KPLBmPP1
      real kTnIPKA,KTnImPKA,kTnIPP2A,KTnImPP2A,PP2Acyt
      real kco,koc,kcop,alphap,PP1cav,PP2Acav,PPcav,ICaLtot,ICaLcav
      real ICaLecav,kICaLPKA,KICaLmPKA,kICaLPP,KICaLmPP,SFICAp,fICaLcav
      real Kapp,Kamp,Kbpp,Kbmp,Kcpp,Kcmp
      real RyRtot,RyRecav,kRyRPKA,kRyRPP,KRyRmPKA,KRyRmPP,f_RyR
      real INatot,INacav,kINaPKA,kINaPP,KINamPKA,KINamPP,Gnap,Vap,Vip
      real alp11p,alp12p,alp13p,bet2p
      real kINaKPKA,kINaKPP,KINaKmPKA,KINaKmPP
      real Gkur2p,PP1ecav,kIKurPKA,kIKurPP,KIKurmPKA,KIKurmPP
      real kIK1PKA,kIK1PP,KIK1mPKA,KIK1mPP,tempa7p,tempa8p,Gk1p
      real kIKtoPKA,kIKtoPP,KIKtomPKA,KIKtomPP,Gktop
      real Icat,Icatc,Icats,k_CaTpv,k_CaTmv,k_CaTpo,k_CaTmo,k_CaTpi
      real k_CaTmi,f_CaT,h_CaT,g_CaT,g_CaTp,fICaTcav,kIcatPKA,kIcatPP
      real KIcatmPKA,KIcatmPP
      real k_CaKf1,k_CaKf2,k_CaKf3,k_CaKb1,k_CaKb2,k_CaKb3,k_CaKfo1
      real k_CaKbo1,k_CaKfo2,k_CaKbo2,g_CaK,ICaK
      real kIcakPKA,kIcakPP,KIcakmPKA,KIcakmPP
      real Y(300),YDOT(300)
      real dito
*
      common/a1/Jrelt,Jtrt,Jxfert,Jleak,Jup,Jtrpn
      common/a2/Acap,Vmyo,Vjsr,Vnsr,Vss,Cm,F,Bi,Bss,Bjsr
      common/a10/Icas,Icab,Inaca,Ipca,Ina,Inab,Iclca,Icasc,Icase
      common/a11/Inak,Ikto,Ik1,Iks,Ikur,Ikr,Ikur1,Ikur2,Ikur3,Ikur4
      common/b1/RCcavf,RCecavf,RCcytf,Gscavabg,Gsecavabg,Gscytabg
      common/b2/cAMPcell,Ccell,ACactcell,PDEactcell
      common/b3/fluxcavecav,fluxcavcyt,fluxecavcav
      common/b4/fluxecavcyt,fluxcytcav,fluxcytecav 
      common/b5/RB1cavtot,RB1ecavtot,RB1cyttot
      common/b6/PDE2cavtot,PDE3cavtot,PDE4cavtot,PDE2ecavtot,
     *  PDE3ecavtot,PDE4ecavtot,PDE2cyttot,PDE3cyttot,PDE4cyttot,
     *  PDE8cavtot,PDE8ecavtot,PDE8cyttot
      common/b7/PDEact2,PDEact3,PDEact4,PDEact8,PKAactcell
      common/b8/PKAcav,PKAecav,PKAcyt,PKIcavf,PKIecavf,PKIcytf
      common/b9/RB1pPKA,RB1pbARK,RB1ptot,L,RB2pPKA,RB2pbARK,RB2ptot
      common/b10/PDE2mem,PDE3mem,PDE4mem,PDE8mem,PDEmem
      common/b11/pp1cytf,fICaLcav
      common/b12/Gicavabg,Giecavabg,RB2cavtot,RB2ecavtot
      common/c1/RB2GicavPKA,LRB2GicavPKA,RB2cavPKAf,Gicavf,
     * acavB2i,bcavB2i,ccavB2i
      common/c2/RB2GiecavPKA,LRB2GiecavPKA,RB2ecavPKAf,Giecavf,
     * aecavB2i,becavB2i,cecavB2i
      common/d1/Icat,Icatc,Icats,ICaK
*
* Cell Geometry Parameters
*
      Acap =  0.5931e-4
      Vcell= 22.1200e-6
      Vmyo = 13.9600e-6
      Vcyt = 13.9600e-6
      Vjsr =  0.1349e-6
      Vnsr =  2.5523e-6
      Vss  =  1.4850e-9*2.0
      Vcav =  0.02*Vcell
      Vecav=  0.04*Vcell
*
* Standard Ionic Concentrations
*
      cKo  =   5400.0
      cNao = 140000.0
      cCao =   1800.0
*
* SR Parameters
*
      Icasmax = 2.0
      t1  = 0.04
      t2  = -0.1/Icasmax
      v1    =  4.5*0.09
      v2    =  1.74e-5*3.0
      v3    =  0.306*2.63
*      Kmup  =  0.5
      Kmup  = (0.38*(1-Y(104))+0.25*Y(104))*0.815
      Ttr   = 20.0
      Txfer =  8.0*1.3
      Kap   =  0.006075
      Kam   =  0.07125
      Kbp   =  0.00405
      Kbm   =  0.965
      Kcp   =  0.009
      Kcm   =  0.0008
      Kapp  =  Kap*5.0
      Kamp  =  Kam*3.0
      Kbpp  =  Kbp*5.0
      Kbmp  =  Kbm*3.0
      Kcpp  =  Kcp*50.0
      Kcmp  =  Kcm*30.0
      f_RyR =  0.001
*
* L-type Ca2+ Channel Parameters
*
      Kpcmax  =  0.11662*2.0
      Kpchalf = 10.0
      Kpcb    =  0.0005*4.0*1.2
      Kpcf    =  2.5*8.0*2.0
      kco     =  1.0
      kcop    =  4.0
      koc     =  1.0
      SFC4I1  =  0.01
      SFC4I2  =  0.002
      SFC4I3  =  1.0
      SFOI1   =  1.0
      SFOI2   =  0.001
      SFI1I3  =  0.001
      SFI2I3  =  1.0
      SFICA   =  0.1127
      SFICAp  =  0.2351
*
* T-type Ca2+ Channel Parameters
*
      k_CaTpv = 9.5*exp(Y(1)/15.0)
      k_CaTmv = 0.008*exp(-1.0*Y(1)/18.0)*4.0
      k_CaTpo = 4.0
      k_CaTmo = 0.025*exp(-1.0*Y(1)/34.0)
      k_CaTpi = 0.070
      k_CaTmi = 0.0014*0.7
      f_CaT   = 0.2
      h_CaT   = 0.5
      g_CaT   = 0.012725
      g_CaTp  = 0.02545
      fICaTcav= 0.68
*
* Buffering Parameters
*
      cLTRPNtot =    70.0
      cHTRPNtot =   140.0
      khtrpnp   =     0.00237
      khtrpnm   =     0.000032
      kltrpnp   =     0.0327
*      kltrpnm   =     0.0196
      kltrpnm   =     0.0196*(1-Y(105))+0.0294*Y(105)
      cCMDNtot  =    50.0
      cCSQNtot  = 15000.0*0.6
      Kmcmdn    =     0.238
      Kmcsqn    =   800.0
*
* Membrane Current Parameters
*
      Cm      =     1.0
      F       =    96.5
      T       =   298.0
      R       =     8.314
      factor  =    R*T/F
      ifactor =    F/(R*T)
      Gna     =    14.4
      Gnap    =    18.0
      knaca   =   275.0*1.0
      Kmna    = 87500.0
      Kmca    =  1380.0
      ksat    =     0.27
      nu      =     0.35
      Inakmax =     4.0*0.5
*      Kmnai   = 21000.0
      Kmnai   = 18800.0*(1-Y(129))+13600.0*Y(129)
      Kmko    =  1500.0
      sigma   = (exp(cNao/67300.0)-1.0)/7.0
      Ipcamax =     0.051*1.5*2.0
      Kmpca   =     0.5
      Gcab    =     0.000284*0.75
      Gnab    =     0.0063*1.0
      dito    =     7.0
*      Gks     =     0.00575
      Gks     =     0.0
      Gkto    =     0.14067
      Gktop   =     0.14067*(1.0-0.387)
      Gkur1   =     0.0
      Gkur2   =     0.05766
      Gkur2p  =     0.07496
      Gkur3   =     0.0428
      Gkur4   =     0.05766
*
* HERG current parameters
*
      kf  = 0.023761
      kb  = 0.036778
      Gkr = 0.078*10.0
*
* Calcium activated chloride current
*
      Gclca     = 10.0
      Poclcamax = 0.2
      Kmclca    = 10.0
      Ecl       = -40.0
*
* Small-conductance calcium-activated K+ current
*
      k_CaKf1    = 0.2*Y(2)*0.6
      k_CaKf2    = 0.16*Y(2)*0.6
      k_CaKf3    = 0.08*Y(2)*0.6
      k_CaKb1    = 0.08
      k_CaKb2    = 0.08
      k_CaKb3    = 0.2
      k_CaKfo1   = 0.16
      k_CaKbo1   = 1.0
      k_CaKfo2   = 1.2
      k_CaKbo2   = 0.1
      g_CaK      = (1.1E-6*Y(1)*Y(1)*Y(1)+1.0E-4*Y(1)*Y(1)
     *           + 0.0041*Y(1)+0.1872)*4.0
*
* Beta1-adrenergic receptor module
*
*      L       =   0.0
      RB1tot  =   0.0103
      fcavB1  =   0.01
      fecavB1 =   0.5
      fcytB1  =   1-fcavB1-fecavB1
      Gstot   =   2.054
      fcavGs  =   0.4
      fecavGs =   0.4
      fcytGs  =   1-fcavGs-fecavGs
      KB1L    =   0.567
      KB1H    =   0.0617
      KB1C    =   2.86
      kPKAp   =   8.1e-7
      kPKAm   =   0.25*kPKAp
      kBARKp  =   0.3*kPKAp
      kBARKm  =   kPKAm
      factGsGi=   0.04
      kact1Gs =   4.9e-3
      kact2Gs =   2.6e-4
      khydGs  =   8.0e-4
      kreasGs =   1.2
*
* Beta2-adrenergic receptor module
*
      pi      =   3.14159026
      RB2tot  =   0.0053
      fcavB2  =   0.99
      fecavB2 =   1-fcavB2
      Gitot   =  10.086
      fcavGi  =   0.99
      fecavGi =   1-fcavGi
      KB2L    =   1.053
      KB2H    =   0.0118
      KB2C    =   5.86
      KB2F    =   0.0189
      KB2A    =  28.79
      kact1Gi =   4.0e-3*0.5
      kact2Gi =   5.0e-5
      khydGi  =   khydGs
      kreasGi =   kreasGs
*
* Adenylyl Cyclase module
*   
      Kmatp  = 340.0
      ATP    = 5000.0
      ACtot  = 0.02622
      fac56ac47 = 0.74
      fcavac56  = 0.0875
      fecavac47 = 0.1648
      Kac56mgsa = 0.0852
      hac56gsa =  1.357
      Vac56gbg  = 1.430 
      Kac56mgsbg= 0.003793
      hac56gsbg = 1.0842  
      AC56bas   = 0.0377
      AF56      = 0.0511335
      Kac47mgsa = 0.05008 
      hac47gsa  = 1.1657  
      Vac47gbg  = 1.35
      Kac47mgsbg= 0.004466
      hac47gsbg = 0.870   
      AC47bas   = 0.04725
      AF47      = 0.009283
      Kac56mgsgi= 0.482
      Kac56mgi  = 0.0465
      hac56gsgi = 0.662
      Vac56gsgi = 0.857
*
* Phosphodiesterase Module
*
      IBMX = 1.0e-20
      hibmxpde2 = 1.0
      Kibmxmpde2 = 29.5
      hibmxpde3 = 1.0
      Kibmxmpde3 = 5.1
      hibmxpde4 = 1.0
      Kibmxmpde4 = 16.2
      kfpdep    = 1.96e-5
      kbpdep    = 1.02e-5
      deltakpde34 = 3.0
      kpde2     = 0.020
      Kmpde2    = 33.0
      kpde3     = 0.0025
      Kmpde3    = 0.44
      kpde4     = 0.0035
      Kmpde4    = 1.4
      kpde8     = 0.0
      Kmpde8    = 0.15
      fpdepart  = 0.2
      rpartpde23= 0.570
      rpartpde34= 0.748
      PDE2tot   = 0.034610
      PDE3tot   = 0.010346
      PDE4tot   = 0.026687
      PDE8tot   = 0.016805
      fcavpde2  = 0.06608
      fecavpde2 = 2*fcavpde2
      fcytpde2  = 1-fcavpde2-fecavpde2
      fcavpde3  = 0.29814
      fecavpde3 =  0.0
      fcytpde3  = 1-fcavpde3-fecavpde3
      fcavpde4  = 0.05366
      fecavpde4 = 2*fcavpde4
      fcytpde4  = 1-fcavpde4-fecavpde4
      fcavpde8  = 0.06608
      fecavpde8 = 2*fcavpde8
      fcytpde8  = 1-fcavpde8-fecavpde8
*
* cAMP-PKA module
*
      PKAtot   = 0.5176
      fcavpka  = 0.08
      fecavpka = 0.20
      fcytpka  = 1-fcavpka-fecavpka
      PKItot   = 2*0.2*PKAtot
      fcavpki  = fcavpka
      fecavpki = fecavpka
      fcytpki  = fcytpka
      kpkaif1   = 0.0056
      kpkai1    = 2.9
      kpkaif2   = kpkaif1
      kpkai2    = 2.9
*      kpkaif3   = 7.99e-4
      kpkaif3   = 0.0026
      kpkai3    = 1.3
      kpkaiif1  = kpkaif1
      kpkaii1   = 2.5
      kpkaiif2  = kpkaif1
      kpkaii2   = 2.5
*      kpkaiif3  = 0.00511
      kpkaiif3  = kpkaif3
      kpkaii3   = kpkai3
      kpkif     = 0.050
      kpki      = 2.6e-4
*
* Protein Phosphatase module
*
      pp1cyttot = 0.2
      inhib1cyttot = 0.08543
      Kinhib1 = 1.0e-3
      kpkainhib1 = 1.08
      Kmpkainhib1 = 1.5
      kpp2ainhib1pp2acyt = 0.00308
      Kmpp2ainhib1 = 0.001
*
* Caveolae
*
      RB1cavtot   = fcavB1*RB1tot*Vcell/Vcav
      RB2cavtot   = fcavB2*RB2tot*Vcell/Vcav
      Gscavabg    = fcavGs*Gstot*Vcell/Vcav-Y(47)-Y(49)
      Gicavabg    = fcavGi*Gitot*Vcell/Vcav-Y(157)-Y(158)
      RB1cavnptot = RB1cavtot-Y(45)-Y(46)
      RB2cavnptot = RB2cavtot-Y(155)-Y(156)
      acavB2i     = (KB2L+L)*(KB2F+L)/KB2L
      bcavB2i     = Gicavabg*(KB2F+L)-Y(155)*(KB2F+L)
     *            + KB2A*KB2F*(1.0+L/KB2L)
      ccavB2i     =-1.0*Y(155)*KB2A*KB2F
      RB2cavPKAf  =(-bcavB2i+sqrt(bcavB2i**2-4.0*acavB2i*ccavB2i))
     *             /(2.0*acavB2i)
      Gicavf      = Gicavabg/(1.0+RB2cavPKAf*(1.0/KB2A+L/(KB2A*KB2F)))
      LRB2cavPKA  = L*RB2cavPKAf/KB2L
      RB2GicavPKA = RB2cavPKAf*Gicavf/KB2A
      LRB2GicavPKA = L*RB2cavPKAf*Gicavf/(KB2A*KB2F)
      acavB2s     = (KB1H+L)*(KB2H+L)
      bcavB2s     = (L+KB1H)*(L+KB2H)*(RB1cavnptot+RB2cavnptot)
     *            + (KB1C*KB1H+L*KB1C*KB1H/KB1L)*(L+KB2H)
     *            + (KB2C*KB2H+L*KB2C*KB2H/KB2L)*(L+KB1H)
     *            - Gscavabg*(L+KB1H)*(L+KB2H)
      ccavB2s     = (L+KB1H)*(KB2C*KB2H+L*KB2C*KB2H/KB2L)
     *            *(RB1cavnptot-Gscavabg)
     *            + (L+KB2H)*(KB1C*KB1H+L*KB1C*KB1H/KB1L)
     *            *(RB2cavnptot-Gscavabg)
     * +(KB1C*KB1H+L*KB1C*KB1H/KB1L)*(KB2C*KB2H+L*KB2C*KB2H/KB2L)
      dcavB2s     = -Gscavabg*(KB1C*KB1H+L*KB1C*KB1H/KB1L)
     *            *(KB2C*KB2H+L*KB2C*KB2H/KB2L)
*
* Cubic equation solution
*
      pcav        = bcavB2s/acavB2s
      qcav        = ccavB2s/acavB2s
      rcav        = dcavB2s/acavB2s
      Acav        = (3.0*qcav-pcav*pcav)/3.0
      Bcav        = (2.0*pcav*pcav*pcav-9.0*pcav*qcav+27.0*rcav)/27.0
      Dcav        = Acav*Acav*Acav/27.0+Bcav*Bcav/4.0
      Mcav        = (-Bcav*0.5+sqrt(Dcav))**(1.0/3.0)
      Ncav        = (-Bcav*0.5-sqrt(Dcav))**(1.0/3.0)
*
      if (Dcav .gt. 0.0) then
      y1cav = Mcav + Ncav
      y2cav = 0.0
      y3cav = 0.0
      end if
      if (Dcav .eq. 0.0) then
      y1cav = Mcav + Ncav
      y2cav = (Mcav + Ncav)/(-2.0)
      y3cav = (Mcav + Ncav)/(-2.0)
      end if
      if (Dcav .lt. 0.0) then
        if (Bcav .gt. 0.0) then
        ficav = acos(-sqrt((Bcav*Bcav/4.0)/(-Acav*Acav*Acav/27.0)))
        else
        ficav = acos(sqrt((Bcav*Bcav/4.0)/(-Acav*Acav*Acav/27.0)))
        end if
      y1cav = 2.0*sqrt(-Acav/3.0)*cos(ficav)
      y2cav = 2.0*sqrt(-Acav/3.0)*cos(ficav+2.0*pi/3.0)
      y3cav = 2.0*sqrt(-Acav/3.0)*cos(ficav+4.0*pi/3.0)
      end if
*
      z1cav = y1cav - pcav/3.0
      z2cav = y2cav - pcav/3.0
      z3cav = y3cav - pcav/3.0
      Gscavf = max(z1cav,z2cav,z3cav)
      RB1cavnpf = RB1cavnptot
     * /(1.0+L/KB1L+Gscavf*(L/(KB1C*KB1H)+1.0/KB1C))
      RB2cavnpf = RB2cavnptot
     * /(1.0+L/KB2L+Gscavf*(L/(KB2C*KB2H)+1.0/KB2C))
      LRB1cavnp   = L*RB1cavnpf/KB1L
      RB1Gscavnp  = RB1cavnpf*Gscavf/KB1C
      LRB1Gscavnp = L*RB1cavnpf*Gscavf/(KB1C*KB1H)
      LRB2cavnp   = L*RB2cavnpf/KB2L
      RB2Gscavnp  = RB2cavnpf*Gscavf/KB2C
      LRB2Gscavnp = L*RB2cavnpf*Gscavf/(KB2C*KB2H)
      AC56cav   = fcavac56*fac56ac47*ACtot*(Vcell/Vcav)   
      kcavAC56  = AF56*(AC56bas+Y(47)**hac56gsa/(Kac56mgsa
     *          +Y(47)**hac56gsa))*(1+Vac56gbg*Y(48)**hac56gsbg
     *          /(Kac56mgsbg+Y(48)**hac56gsbg))
     *          *(1.0-(1.0-Vac56gsgi*Y(47)**hac56gsgi/(Kac56mgsgi
     *          +Y(47)**hac56gsgi))*Y(157)/(Kac56mgi+Y(157)))
      PDE2cavtot = (1.0-IBMX**hibmxpde2/(Kibmxmpde2+IBMX**hibmxpde2))
     *          *fcavpde2*PDE2tot*(Vcell/Vcav)
      PDE3cavtot = (1.0-IBMX**hibmxpde3/(Kibmxmpde3+IBMX**hibmxpde3))
     *          *fcavpde3*PDE3tot*(Vcell/Vcav)
      PDE4cavtot = (1.0-IBMX**hibmxpde4/(Kibmxmpde4+IBMX**hibmxpde4))
     *          *fcavpde4*PDE4tot*(Vcell/Vcav)
      PDE8cavtot = fcavpde8*PDE8tot*(Vcell/Vcav)
      PKAcav     = fcavpka*PKAtot*(Vcell/Vcav)
      RCcavf     = 2.0*PKAcav-Y(80)-Y(81)-Y(82) 
      PKIcavf    = fcavpki*PKItot*(Vcell/Vcav)-Y(84)
      kpkaiib1   = kpkaiif1*kpkaii1
      kpkaiib2   = kpkaiif2*kpkaii2
      kpkaiib3   = kpkaiif3/kpkaii3
      kpkib      = kpkif*kpki
*
* Extracaveolae
*
      RB1ecavtot   = fecavB1*RB1tot*Vcell/Vecav
      RB2ecavtot   = fecavB2*RB2tot*Vcell/Vecav
      Gsecavabg    = fecavGs*Gstot*Vcell/Vecav-Y(52)-Y(54)
      Giecavabg    = fecavGi*Gitot*Vcell/Vecav-Y(161)-Y(162)
      RB1ecavnptot = RB1ecavtot-Y(50)-Y(51)
      RB2ecavnptot = RB2ecavtot-Y(159)-Y(160)
      aecavB2i     = (KB2L+L)*(KB2F+L)/KB2L
      becavB2i     = Giecavabg*(KB2F+L)-Y(159)*(KB2F+L)
     *             + KB2A*KB2F*(1.0+L/KB2L)
      cecavB2i     =-1.0*Y(159)*KB2A*KB2F
      RB2ecavPKAf  =(-becavB2i+sqrt(becavB2i**2
     *             -4.0*aecavB2i*cecavB2i))/(2.0*aecavB2i)
      Giecavf      = Giecavabg/(1.0+RB2ecavPKAf
     *             *(1.0/KB2A+L/(KB2A*KB2F)))
      LRB2ecavPKA  = L*RB2ecavPKAf/KB2L
      RB2GiecavPKA = RB2ecavPKAf*Giecavf/KB2A
      LRB2GiecavPKA = L*RB2ecavPKAf*Giecavf/(KB2A*KB2F)
      aecavB2s     = (KB1H+L)*(KB2H+L)
      becavB2s     = (L+KB1H)*(L+KB2H)*(RB1ecavnptot+RB2ecavnptot)
     *             + (KB1C*KB1H+L*KB1C*KB1H/KB1L)*(L+KB2H)
     *             + (KB2C*KB2H+L*KB2C*KB2H/KB2L)*(L+KB1H)
     *             - Gsecavabg*(L+KB1H)*(L+KB2H)
      cecavB2s     = (L+KB1H)*(KB2C*KB2H+L*KB2C*KB2H/KB2L)
     *             *(RB1ecavnptot-Gsecavabg)
     *             + (L+KB2H)*(KB1C*KB1H+L*KB1C*KB1H/KB1L)
     *             *(RB2ecavnptot-Gsecavabg)
     * +(KB1C*KB1H+L*KB1C*KB1H/KB1L)*(KB2C*KB2H+L*KB2C*KB2H/KB2L)
      decavB2s     = -Gsecavabg*(KB1C*KB1H+L*KB1C*KB1H/KB1L)
     *             *(KB2C*KB2H+L*KB2C*KB2H/KB2L)
*
* Cubic equation solution
*
      pecav        = becavB2s/aecavB2s
      qecav        = cecavB2s/aecavB2s
      recav        = decavB2s/aecavB2s
      Aecav        = (3.0*qecav-pecav*pecav)/3.0
      Becav        = (2.0*pecav*pecav*pecav-9.0*pecav*qecav
     *             +27.0*recav)/27.0
      Decav        = Aecav*Aecav*Aecav/27.0+Becav*Becav/4.0
      Mecav        = (-Becav*0.5+sqrt(Decav))**(1.0/3.0)
      Necav        = (-Becav*0.5-sqrt(Decav))**(1.0/3.0)
*      
      if (Decav .gt. 0.0) then
      y1ecav = Mecav + Necav
      y2ecav = 0.0
      y3ecav = 0.0
      end if
      if (Decav .eq. 0.0) then
      y1ecav = Mecav + Necav
      y2ecav = (Mecav + Necav)/(-2.0)
      y3ecav = (Mecav + Necav)/(-2.0)
      end if
      if (Decav .lt. 0.0) then
        if (Becav .gt. 0.0) then
        fiecav = acos(-sqrt((Becav*Becav/4.0)
     *         /(-Aecav*Aecav*Aecav/27.0)))
        else
        fiecav = acos(sqrt((Becav*Becav/4.0)
     *         /(-Aecav*Aecav*Aecav/27.0)))
        end if
      y1ecav = 2.0*sqrt(-Aecav/3.0)*cos(fiecav)
      y2ecav = 2.0*sqrt(-Aecav/3.0)*cos(fiecav+2.0*pi/3.0)
      y3ecav = 2.0*sqrt(-Aecav/3.0)*cos(fiecav+4.0*pi/3.0)
      end if
*
      z1ecav = y1ecav - pecav/3.0
      z2ecav = y2ecav - pecav/3.0
      z3ecav = y3ecav - pecav/3.0
      Gsecavf = max(z1ecav,z2ecav,z3ecav)
      RB1ecavnpf = RB1ecavnptot
     * /(1.0+L/KB1L+Gsecavf*(L/(KB1C*KB1H)+1.0/KB1C))
      RB2ecavnpf = RB2ecavnptot
     * /(1.0+L/KB2L+Gsecavf*(L/(KB2C*KB2H)+1.0/KB2C))
      LRB1ecavnp   = L*RB1ecavnpf/KB1L
      RB1Gsecavnp  = RB1ecavnpf*Gsecavf/KB1C
      LRB1Gsecavnp = L*RB1ecavnpf*Gsecavf/(KB1C*KB1H)
      LRB2ecavnp   = L*RB2ecavnpf/KB2L
      RB2Gsecavnp  = RB2ecavnpf*Gsecavf/KB2C
      LRB2Gsecavnp = L*RB2ecavnpf*Gsecavf/(KB2C*KB2H)
      AC47ecav     = fecavac47*(1.0-fac56ac47)*ACtot*Vcell/Vecav
      kecavAC47    = AF47*(AC47bas+Y(52)**hac47gsa/(Kac47mgsa
     *             +Y(52)**hac47gsa))*(1.0+Vac47gbg*Y(53)**hac47gsbg
     *             /(Kac47mgsbg+Y(53)**hac47gsbg))
      PDE2ecavtot=(1.0-(IBMX**hibmxpde2/(Kibmxmpde2+IBMX**hibmxpde2)))
     *          *fecavpde2*PDE2tot*(Vcell/Vecav)
      PDE3ecavtot=(1.0-(IBMX**hibmxpde3/(Kibmxmpde3+IBMX**hibmxpde3)))
     *          *fecavpde3*PDE3tot*(Vcell/Vecav)
      PDE4ecavtot=(1.0-(IBMX**hibmxpde4/(Kibmxmpde4+IBMX**hibmxpde4)))
     *          *fecavpde4*PDE4tot*(Vcell/Vecav)
      PDE8ecavtot = fecavpde8*PDE8tot*(Vcell/Vecav)
      PKAecav     = fecavpka*PKAtot*(Vcell/Vecav)
      RCecavf     = 2.0*PKAecav-Y(86)-Y(87)-Y(88)
      PKIecavf    = fecavpki*PKItot*(Vcell/Vecav)-Y(90)
      kpkaiib1   = kpkaiif1*kpkaii1
      kpkaiib2   = kpkaiif2*kpkaii2
      kpkaiib3   = kpkaiif3/kpkaii3
      kpkib      = kpkif*kpki
*
* Cytosol
*
      RB1cyttot   = fcytB1*RB1tot*Vcell/Vcyt
      Gscytabg    = fcytGs*Gstot*Vcell/Vcyt-Y(57)-Y(59)
      RB1cytnptot = RB1cyttot-Y(55)-Y(56)
      acytB1      = (KB1L+L)*(KB1H+L)/KB1L
      bcytB1      = Gscytabg*(KB1H+L)-RB1cytnptot*(KB1H+L)
     *             +KB1C*KB1H*(1.0+L/KB1L)
      ccytB1      =-RB1cytnptot*KB1C*KB1H
      RB1cytnpf   =(-bcytB1+sqrt(bcytB1**2-4.0*acytB1*ccytB1))
     *             /(2.0*acytB1)
      Gscytf      = Gscytabg/(1.0+RB1cytnpf*(1.0/KB1C+L/(KB1C*KB1H)))
      LRB1cytnp   = L*RB1cytnpf/KB1L
      RB1Gscytnp  = RB1cytnpf*Gscytf/KB1C
      LRB1Gscytnp = L*RB1cytnpf*Gscytf/(KB1C*KB1H)
      AC56cyt     = (1.0-fcavac56)*fac56ac47*ACtot*(Vcell/Vcyt)
      kcytAC56    = AF56*(AC56bas+Y(57)**hac56gsa/(Kac56mgsa
     *            +Y(57)**hac56gsa))*(1.0+Vac56gbg*Y(58)**hac56gsbg
     *            /(Kac56mgsbg+Y(58)**hac56gsbg))
      AC47cyt     = (1.0-fecavac47)*(1.0-fac56ac47)*ACtot*Vcell/Vcyt
      kcytAC47    = AF47*(AC47bas+Y(57)**hac47gsa/(Kac47mgsa
     *            +Y(57)**hac47gsa))*(1.0+Vac47gbg*Y(58)**hac47gsbg
     *            /(Kac47mgsbg+Y(58)**hac47gsbg))
      PDE2cyttot =(1.0-(IBMX**hibmxpde2/(Kibmxmpde2+IBMX**hibmxpde2)))
     *            *fcytpde2*PDE2tot*(Vcell/Vcyt)
      PDE3cyttot =(1.0-(IBMX**hibmxpde3/(Kibmxmpde3+IBMX**hibmxpde3)))
     *            *fcytpde3*PDE3tot*(Vcell/Vcyt)
      PDE4cyttot =(1.0-(IBMX**hibmxpde4/(Kibmxmpde4+IBMX**hibmxpde4)))
     *            *fcytpde4*PDE4tot*(Vcell/Vcyt)
      PDE8cyttot = fcytpde8*PDE8tot*(Vcell/Vcyt)
      PKAcyt      = fcytpka*PKAtot*(Vcell/Vcyt)
      RCcytf      = 2.0*PKAcyt-Y(92)-Y(93)-Y(94)
      PKIcytf     = fcytpki*PKItot*(Vcell/Vcyt)-Y(96)
      kpkaib1     = kpkaif1*kpkai1
      kpkaib2     = kpkaif2*kpkai2
      kpkaib3     = kpkaif3/kpkai3
      kpkib       = kpkif*kpki
*
* Protein phosphatase module
*
      inhib1cytf = inhib1cyttot-Y(100)
      ainh1  = 1.0
      binh1  = Kinh1+pp1cyttot-Y(100)
      cinh1  = -Y(100)*Kinhib1
      inhib1cytp = -binh1/2.0+sqrt(binh1**2-4.0*ainh1*cinh1)/2.0
      pp1cytf = pp1cyttot*Kinhib1/(Kinhib1+inhib1cytp)
*
* cAMP fluxes
*        
      jcavecav = 5.00e-12
      jcavcyt  = 7.500e-11
      jecavcyt = 9.00e-12
*
* PLB module
*
      kPLBPKA  = 1.08917e-4
      KPLBmPKA = 4.90970e-4
      kPLBPP1  = 4.41956e-5
      KPLBmPP1 = 1.69376e-2
*
* TnI module
*
      PP2Acyt   = 0.0607843
      kTnIPKA   = 2.47254e-5
      KTnImPKA  = 2.71430e-5
      kTnIPP2A  = 8.65898e-5
      KTnImPP2A = 8.01420e-1
*
* ICaL module
*
      PP1cav    = 0.1
      PP2Acav   = 0.1
      PPcav     = PP1cav + PP2acav
      ICaLtot   = 0.0273
      fICaLcav  = 0.2
      ICaLcav   = fICaLcav*ICaLtot*Vcell/Vcav
      ICaLecav  = (1.0-fICaLcav)*ICaLtot*Vcell/Vecav
      kICaLPKA  = 2.0e-5*0.87
      kICaLPP   = 1.55e-7*1.5
      KICaLmPKA = 0.5
      KICaLmPP  = 0.2
*
* ICaT module
*
      kIcatPKA  = 1.74e-5*0.5*20.0
      kIcatPP   = 2.33e-7*200.0
      KIcatmPKA = 5.0
      KIcatmPP  = 0.1
*
* ICaK module
*
      kIcakPKA  = 1.74e-5*0.085
      kIcakPP   = 2.33e-7*1.25
      KIcakmPKA = 2.5
      KIcakmPP  = 0.01
*
* RyR module
*
      RyRtot   = 0.1993
      RyRecav  = RyRtot*Vcell/Vecav
      kRyRPKA  = 3.85e-5*1.5
      kRyRPP   = 1.925e-4*1.5
      KRyRmPKA = 0.5
      KRyRmPP  = 0.05
*
* INa module
*
      INatot   = 0.0
      INacav   = INatot*Vcell/Vcav
      kINaPKA  = 6.8400e-6
      kINaPP   = 1.9804e-5
      KINamPKA = 5.49415e-3
      KINamPP  = 0.393025
*
* INaK module
*
      kINaKPKA  = 3.053e-6
      kINaKPP   = 1.8491e-5
      KINaKmPKA = 1.1001e-3
      KINaKmPP  = 5.7392
*
* IKur module
*
      PP1ecav   = 0.1
      kIKurPKA  = 6.9537e-6
      kIKurPP   = 3.1700e-5
      KIKurmPKA = 0.138115
      KIKurmPP  = 0.23310
*
* IK1 module
*
      kIK1PKA  = 1.9954e-5
      kIK1PP   = 9.0968e-5
      KIK1mPKA = 0.027623
      KIK1mPP  = 0.023310
*
* IKto module
*
      kIKtoPKA  = 4.38983e-5
      kIKtoPP   = 9.09678e-5
      KIKtomPKA = 0.27623
      KIKtomPP  = 0.23310
*
      cAMPcell = (Y(97)*Vcav+Y(98)*Vecav+Y(99)*Vcyt)/Vcell
      Ccell    = (Y(83)*Vcav+Y(89)*Vecav+Y(95)*Vcyt)/Vcell
      ACactcell= (YDOT(60)*Vcav+YDOT(61)*Vecav
     *      +(YDOT(62)+YDOT(63))*Vcyt)/Vcell
      PDEact2 = (YDOT(66)*Vcav+YDOT(71)*Vecav+YDOT(76)*Vcyt)/Vcell
      PDEact3 = (YDOT(67)*Vcav+YDOT(72)*Vecav+YDOT(77)*Vcyt)/Vcell
      PDEact4 = (YDOT(68)*Vcav+YDOT(73)*Vecav+YDOT(78)*Vcyt)/Vcell
      PDEact8 = (YDOT(101)*Vcav+YDOT(102)*Vecav+YDOT(103)*Vcyt)/Vcell
      PDE2mem = (YDOT(66)*Vcav+YDOT(71)*Vecav)/Vcell
      PDE3mem = (YDOT(67)*Vcav+YDOT(72)*Vecav)/Vcell
      PDE4mem = (YDOT(68)*Vcav+YDOT(73)*Vecav)/Vcell
      PDE8mem = (YDOT(101)*Vcav+YDOT(102)*Vecav)/Vcell
      PDEactcell = PDEact2+PDEact3+PDEact4+PDEact8
      PDEmem = PDE2mem+PDE3mem+PDE4mem+PDE8mem
      PKAactcell = (YDOT(83)*Vcav+YDOT(89)*Vecav+YDOT(95)*Vcyt)/Vcell
      RB1pPKA = (Y(45)*Vcav+Y(50)*Vecav+Y(55)*Vcyt)/Vcell
      RB1pbARK = (Y(46)*Vcav+Y(51)*Vecav+Y(56)*Vcyt)/Vcell
      RB1ptot = RB1pPKA + RB1pbARK
      RB2pPKA = (Y(155)*Vcav+Y(159)*Vecav)/Vcell
      RB2pbARK = (Y(156)*Vcav+Y(160)*Vecav)/Vcell
      RB2ptot = RB2pPKA + RB2pbARK
*
************************************************************************
*            Calcium fluxes
************************************************************************
*      Pryr     = exp((Y(1)-5.0)*(Y(1)-5.0)/(-648.))
      Jrelt    = v1*(Y(18)+Y(19)+Y(118)+Y(119))*(Y(6)-Y(7))*Y(44)
      Jtrt     = (Y(3)-Y(6))/Ttr
      Jxfert   = (Y(7)-Y(2))/Txfer
      Jleak    = v2*(Y(3)-Y(2))
      Jup      = v3*Y(2)*Y(2)/(Kmup*Kmup+Y(2)*Y(2))
      temp12    = khtrpnp*Y(2)*(cHTRPNtot-Y(5))-khtrpnm*Y(5)
      temp13    = kltrpnp*Y(2)*(cLTRPNtot-Y(4))-kltrpnm*Y(4)
      Jtrpn    = temp12 + temp13
************************************************************************
*             Ionic currents
************************************************************************
*
* Y(1)   : Membrane potential (mV)
*
      V = Y(1)
*
* Icab   : Calcium background current
*
      Ecan = 0.5*factor*alog(cCao/Y(2))
      Icab = Gcab*(Y(1)-Ecan)
*
* Ipca   : Calcium pump current
*
      Ipca = Ipcamax*Y(2)*Y(2)/(Kmpca*Kmpca+Y(2)*Y(2))
*
* Inaca  : Na-Ca exchange current
*
      tempa1 = knaca/(Kmna*Kmna*Kmna+cNao*cNao*cNao)
      tempa2 = 1.0/(Kmca+cCao)
      tempa3 = 1.0/(1.0+ksat*exp((nu-1.0)*Y(1)*ifactor))
          tempa41 = exp(nu*Y(1)*ifactor)*Y(23)*Y(23)*Y(23)*cCao
          tempa42 = exp((nu-1.0)*Y(1)*ifactor)*cNao*cNao*cNao*Y(2)
      tempa4 = tempa41-2.0*tempa42
      Inaca = tempa1*tempa2*tempa3*tempa4
*
* Ica    : L-type calcium current
*
      Icasc = fICaLcav*(SFICA*Y(8)+SFICAp*Y(107))*(Y(1)-52.0)
      Icase = (1.0-fICaLcav)*(SFICA*Y(137)+SFICAp*Y(146))*(Y(1)-52.0)
      Icas  = Icasc + Icase
*
* Icat   : T-type calcium current
*
      Icatc = fICaTcav*(g_CaT*Y(168)*(1-Y(175))+g_CaTp*Y(168)*Y(175))
     *      *(Y(1)-35.0)
      Icats = (1.0-fICaTcav)*(g_CaT*Y(168)*(1-Y(176))
     *      +g_CaTp*Y(168)*Y(176))*(Y(1)-35.0)
      Icat  = Icatc + Icats
*
* Ina    : Na fast current (Luo and Rudy, 1994)
*
      ENa = factor*alog((0.9*cNao+0.1*cKo)/(0.9*Y(23)+0.1*Y(24)))
      Ina = (Gna*Y(37)+Gnap*Y(123))*(Y(1)-ENa)
*
* Inab   : Na background current
*
      Inab = Gnab*(Y(1)-ENa)
*
* Inak    : Na-K exchange current
*
      tempa5 = 1.0/(1.0+(Kmnai/Y(23))**3.0)
      tempa6 = cKo/(cKo+Kmko)
      tempa11= 1.0+0.1245*exp(-0.1*Y(1)*ifactor)
     *         +0.0365*sigma*exp(-1.0*Y(1)*ifactor)
      Inak = Inakmax*tempa5*tempa6/tempa11
*
* Ikto    : transient outward current (Liu and Rasmusson, 1997)
*
      EK = factor*alog(cKo/Y(24))
      Ikto = Gkto*Y(25)*Y(25)*Y(25)*Y(26)*(Y(1)-EK)*(1-Y(134))
     *     + Gktop*Y(135)*Y(135)*Y(135)*Y(136)*(Y(1)-EK)*Y(134)
*
* Ik1    : Time independent K+ current (Rasmusson et al. 1990)
*
*      tempa7 = 0.3397*cKo/(cKo+210.0)
*      tempa8 = 1.0+exp(0.0448*(Y(1)-EK))
*      tempa7p = 0.11394*cKo/(cKo+210.0)
*      tempa8p = 1.0+exp(0.0298*(Y(1)-EK))
*      Ik1 = (tempa7*Y(133)/tempa8
*     *    + tempa7*(1-Y(133))/tempa8)*(Y(1)-EK)
*
      tempa7 = 1.02/(1.0+exp(0.2385*(Y(1)-EK-59.215)))
      tempa8 = (0.8*exp(0.08032*(Y(1)-EK+5.476))
     * +exp(0.06175*(Y(1)-EK-594.31)))
     * /(1.0+exp(-0.5143*(Y(1)-EK+4.753)))
      Gk1  = 1.7*0.27*sqrt(cKo/5400.0)
      Gk1p = 1.7*0.27*sqrt(cKo/5400.0)
      Ik1  = (Gk1*Y(133)*tempa7/(tempa7+tempa8)
     *     + Gk1p*(1-Y(133))*tempa7/(tempa7+tempa8))*(Y(1)-EK)
*
* Iks    : Delayed rectifier K+ current (Rasmusson et al. 1990)
*
      Iks = Gks*Y(27)*Y(27)*(Y(1)-EK)
*
* Ikur   : Ultra-rapidly activating delayed rectifier Ikur
*                  (Zhou et al., 1998)
*
      Ikur1 = Gkur1*Y(28)*Y(29)*(Y(1)-EK)
      Ikur2 = (Gkur2*Y(35)*Y(36)*Y(132)
     *      +  Gkur2p*Y(130)*Y(131)*(1-Y(132)))*(Y(1)-EK)
      Ikur3 = Gkur3*Y(43)*(Y(1)-EK)
      Ikur4 = Gkur4*Y(184)*Y(185)*(Y(1)-EK)
      Ikur = Ikur1+Ikur2+Ikur3+Ikur4
*
* Ikr    : HERG current (Wang et al., 1997)
*
      EKr = factor*alog((0.98*cKo+0.02*cNao)/(0.98*Y(24)+0.02*Y(23)))
      Ikr = Gkr*Y(33)*(Y(1)-EKr)
*
* Iclca  : calcium-activated chloride current (Xu et al., 2002)
*
      tempa9  = Poclcamax/(1+exp((46.7-Y(1))/7.8))
      tempa10 = Y(2)/(Y(2)+Kmclca)
      Iclca = Gclca*tempa9*tempa10*(Y(1)-Ecl)
*
* ICaK : small-conductance calcium-activated K+ current
*
      ICaK = g_CaK*(Y(181)+Y(182))*(1-Y(183))
*
************************************************************************
*
* Y(1)  membrane potential
*
      sum_i = Icas+Icab+Inaca+Ipca+Ina+Inab+Iclca
     *       +Inak+Ikto+Ik1+Iks+Ikur+Ikr-Istim+Icat+ICaK
      YDOT(1) = -sum_i/Cm
*
* Y(2)  intracellular calcium Cai
*
      Bi = 1.0/(1+cCMDNtot*Kmcmdn/((Kmcmdn+Y(2))*(Kmcmdn+Y(2))))
      temp1 = 0.5*(Icab-2.0*Inaca+Ipca+Icasc+Icat)*Acap/(Vmyo*F)
      YDOT(2) = Bi*(Jleak + Jxfert - Jup - Jtrpn - temp1)
*
* Y(3)  network SR calcium Cansr
*
      YDOT(3) = (Jup-Jleak)*Vmyo/Vnsr - Jtrt*Vjsr/Vnsr
*
* Y(4)  cLTRPNca
*
      YDOT(4) = kltrpnp*Y(2)*(cLTRPNtot-Y(4)) - kltrpnm*Y(4)
*
* Y(5)  cHTRPNca
*
      YDOT(5) = khtrpnp*Y(2)*(cHTRPNtot-Y(5)) - khtrpnm*Y(5)
*
*     Partial contributions
*
      alpha = 0.4 * exp((Y(1)+15.0)/15.0)
      alphap = 0.4 * exp((Y(1)+15.0+20.0)/15.0)
      beta  = 0.13 * exp(-1.0*(Y(1)+15.0)/18.0)
*
      gammac = Kpcmax*Y(2)/(Kpchalf+Y(2))
*
* Y(6) cCajsr junction SR calcium concentartion
*
      temp2 = (Kmcsqn+Y(6))*(Kmcsqn+Y(6))
      Bjsr  = 1.0/(1.0+cCSQNtot*Kmcsqn/temp2)
      YDOT(6) = Bjsr*(Jtrt-Jrelt)
*
* Y(7) cCass  subspace calcium concentration
*
      temp3 = (Kmcmdn+Y(7))*(Kmcmdn+Y(7))
      Bss   = 1.0/(1.0+cCMDNtot*Kmcmdn/temp3)
      temp4 = Jrelt*Vjsr/Vss-Jxfert*Vmyo/Vss
      YDOT(7) = Bss*(temp4-Icase*Acap/(2.0*Vss*F))
*
*               ICaLcav
*
* Y(8)     O  Ca channel variable
*
      YDOT(8) =
     *   kco*Y(106) - koc*Y(8)
     * + SFOI1*Kpcb*Y(13) - SFOI1*gammac*Y(8)
     * + SFOI2*alpha*Y(14) - SFOI2*Kpcf*Y(8)
     * - kICaLPKA*Y(83)*Y(8)/(KICaLmPKA+ICaLcav*Y(8))
     * + (alpha/alphap)*kICaLPP*PPcav*Y(107)
     * /(KICaLmPP+ICaLcav*Y(107))
*
* Y(9)    C1  Ca channel variable
*
*      YDOT(9) = beta*Y(10)-4.0*alpha*Y(9)
*
* Y(10)   C2  Ca channel variable
*
      YDOT(10) =
     *    4.0*alpha*Y(9) - beta*Y(10)
     *  + 2.0*beta*Y(11) - 3.0*alpha*Y(10)
     * - kICaLPKA*Y(83)*Y(10)/(KICaLmPKA+ICaLcav*Y(10))
     * + (alphap*alphap*kcop/(alpha*alpha*kco))*kICaLPP*PPcav*Y(109)
     * /(KICaLmPP+ICaLcav*Y(109))
*
* Y(11)   C3  Ca channel variable
*
      YDOT(11) =
     *    3.0*alpha*Y(10) - 2.0*beta*Y(11)
     *  + 3.0*beta*Y(12) - 2.0*alpha*Y(11)
     * - kICaLPKA*Y(83)*Y(11)/(KICaLmPKA+ICaLcav*Y(11))
     * + (alphap*kcop/(alpha*kco))*kICaLPP*PPcav*Y(110)
     * /(KICaLmPP+ICaLcav*Y(110))
*
* Y(12)   C4  Ca channel variable
*
      YDOT(12) =
     *    2.0*alpha*Y(11) - 3.0*beta*Y(12)
     *  + 4.0*beta*Y(106) - alpha*Y(12)
     *  + 4.0*SFC4I1*Kpcb*beta*Y(13) 
     *  - SFC4I1*alpha*gammac*(kco/koc)*Y(12)
     *  + 4.0*SFC4I2*beta*Y(14) - SFC4I2*Kpcf*(kco/koc)*Y(12)
     *  + 4.0*SFC4I3*beta*Kpcb*Y(15) 
     *  - SFC4I3*gammac*Kpcf*(kco/koc)*Y(12)
     * - kICaLPKA*Y(83)*Y(12)/(KICaLmPKA+ICaLcav*Y(12))
     * + (kcop/kco)*kICaLPP*PPcav*Y(111)
     * /(KICaLmPP+ICaLcav*Y(111))
*
* Y(106) CP Ca channel variable
*
      YDOT(106) =
     *    alpha*Y(12) - 4.0*beta*Y(106)
     *  + koc*Y(8) - kco*Y(106)
     * - kICaLPKA*Y(83)*Y(106)/(KICaLmPKA+ICaLcav*Y(106))
     * + (alpha*kcop/(alphap*kco))*kICaLPP*PPcav*Y(115)
     * /(KICaLmPP+ICaLcav*Y(115))
*
* Y(13)   I1  Ca channel variable
*
      YDOT(13) =
     *    SFOI1*gammac*Y(8) - SFOI1*Kpcb*Y(13)
     *  + SFI1I3*alpha*Y(15) - SFI1I3*Kpcf*Y(13)
     *  + SFC4I1*alpha*gammac*(kco/koc)*Y(12) 
     *  - 4.0*SFC4I1*beta*Kpcb*Y(13)
     * - kICaLPKA*Y(83)*Y(13)/(KICaLmPKA+ICaLcav*Y(13))
     * + (alpha/alphap)*kICaLPP*PPcav*Y(112)
     * /(KICaLmPP+ICaLcav*Y(112))
*
* Y(14)   I2  Ca channel variable
*
      YDOT(14) =
     *    SFOI2*Kpcf*Y(8) - SFOI2*alpha*Y(14)
     *  + SFI2I3*Kpcb*Y(15) - SFI2I3*gammac*Y(14)
     *  + SFC4I2*Kpcf*(kco/koc)*Y(12) - 4.0*SFC4I2*beta*Y(14)
     * - kICaLPKA*Y(83)*Y(14)/(KICaLmPKA+ICaLcav*Y(14))
     * + kICaLPP*PPcav*Y(113)/(KICaLmPP+ICaLcav*Y(113))
*
* Y(15)   I3  Ca channel variable
*
      YDOT(15) =
     *    SFI1I3*Kpcf*Y(13) - SFI1I3*alpha*Y(15)
     *  + SFI2I3*gammac*Y(14) - SFI2I3*Kpcb*Y(15)
     *  + SFC4I3*gammac*Kpcf*(kco/koc)*Y(12) 
     *  - 4.0*SFC4I3*beta*Kpcb*Y(15)
     * - kICaLPKA*Y(83)*Y(15)/(KICaLmPKA+ICaLcav*Y(15))
     * + kICaLPP*PPcav*Y(114)/(KICaLmPP+ICaLcav*Y(114))
*
* Y(107) OP Ca channel variable
*
      YDOT(107) =
     *   kcop*Y(115) - koc*Y(107)
     * + SFOI1*Kpcb*Y(112) - SFOI1*gammac*Y(107)
     * + SFOI2*alphap*Y(113) - SFOI2*Kpcf*Y(107)
     * + kICaLPKA*Y(83)*Y(8)/(KICaLmPKA+ICaLcav*Y(8))  
     * - (alpha/alphap)*kICaLPP*PPcav*Y(107)
     * /(KICaLmPP+ICaLcav*Y(107))
*
* Y(108) C1P Ca channel variable
*
      YDOT(108) = beta*Y(109)-4.0*alphap*Y(108)
     * + kICaLPKA*Y(83)*Y(9)/(KICaLmPKA+ICaLcav*Y(9))
     * - (alphap*alphap*alphap*kcop/(alpha*alpha*alpha*kco))
     * *kICaLPP*PPcav*Y(108)/(KICaLmPP+ICaLcav*Y(108))
*
* Y(109) C2P Ca channel variable
*
      YDOT(109) =
     *    4.0*alphap*Y(108) - beta*Y(109)
     *  + 2.0*beta*Y(110) - 3.0*alphap*Y(109)
     * + kICaLPKA*Y(83)*Y(10)/(KICaLmPKA+ICaLcav*Y(10))
     * - (alphap*alphap*kcop/(alpha*alpha*kco))*kICaLPP*PPcav*Y(109)
     * /(KICaLmPP+ICaLcav*Y(109))
*
* Y(110) C3P Ca channel variable
*
      YDOT(110) =
     *    3.0*alphap*Y(109) - 2.0*beta*Y(110)
     *  + 3.0*beta*Y(111) - 2.0*alphap*Y(110)
     * + kICaLPKA*Y(83)*Y(11)/(KICaLmPKA+ICaLcav*Y(11))
     * - (alphap*kcop/(alpha*kco))*kICaLPP*PPcav*Y(110)
     * /(KICaLmPP+ICaLcav*Y(110))
*
* Y(111) C4P Ca channel variable
*
      YDOT(111) =
     *    2.0*alphap*Y(110) - 3.0*beta*Y(111)
     *  + 4.0*beta*Y(115) - alphap*Y(111)
     *  + 4.0*SFC4I1*Kpcb*beta*Y(112)
     *  - SFC4I1*alphap*gammac*(kcop/koc)*Y(111)
     *  + 4.0*SFC4I2*beta*Y(113) - SFC4I2*Kpcf*(kcop/koc)*Y(111)
     *  + 4.0*SFC4I3*beta*Kpcb*Y(114)
     *  - SFC4I3*gammac*Kpcf*(kcop/koc)*Y(111)
     * + kICaLPKA*Y(83)*Y(12)/(KICaLmPKA+ICaLcav*Y(12))
     * - (kcop/kco)*kICaLPP*PPcav*Y(111)
     * /(KICaLmPP+ICaLcav*Y(111))
*
* Y(115) CPP Ca channel variable
*
      YDOT(115) =
     *    alphap*Y(111) - 4.0*beta*Y(115)
     *  + koc*Y(107) - kcop*Y(115)
     * + kICaLPKA*Y(83)*Y(106)/(KICaLmPKA+ICaLcav*Y(106))
     * - (alpha*kcop/(alphap*kco))*kICaLPP*PPcav*Y(115)
     * /(KICaLmPP+ICaLcav*Y(115))
*
* Y(112) I1P Ca channel variable
*
      YDOT(112) =
     *    SFOI1*gammac*Y(107) - SFOI1*Kpcb*Y(112)
     *  + SFI1I3*alphap*Y(114) - SFI1I3*Kpcf*Y(112)
     *  + SFC4I1*alphap*gammac*(kcop/koc)*Y(111)
     *  - 4.0*SFC4I1*beta*Kpcb*Y(112)
     * + kICaLPKA*Y(83)*Y(13)/(KICaLmPKA+ICaLcav*Y(13))
     * - (alpha/alphap)*kICaLPP*PPcav*Y(112)
     * /(KICaLmPP+ICaLcav*Y(112))
*
* Y(113) I2P Ca channel variable
*
      YDOT(113) =
     *    SFOI2*Kpcf*Y(107) - SFOI2*alphap*Y(113)
     *  + SFI2I3*Kpcb*Y(114) - SFI2I3*gammac*Y(113)
     *  + SFC4I2*Kpcf*(kcop/koc)*Y(111) - 4.0*SFC4I2*beta*Y(113)
     * + kICaLPKA*Y(83)*Y(14)/(KICaLmPKA+ICaLcav*Y(14))
     * - kICaLPP*PPcav*Y(113)/(KICaLmPP+ICaLcav*Y(113))
*
* Y(114) I3P Ca channel variable
*
      YDOT(114) =
     *    SFI1I3*Kpcf*Y(112) - SFI1I3*alphap*Y(114)
     *  + SFI2I3*gammac*Y(113) - SFI2I3*Kpcb*Y(114)
     *  + SFC4I3*gammac*Kpcf*(kcop/koc)*Y(111)
     *  - 4.0*SFC4I3*beta*Kpcb*Y(114)
     * + kICaLPKA*Y(83)*Y(15)/(KICaLmPKA+ICaLcav*Y(15))
     * - kICaLPP*PPcav*Y(114)/(KICaLmPP+ICaLcav*Y(114))
*
*     Partial contributions
*
      gammae = Kpcmax*Y(7)/(Kpchalf+Y(7))
*
*                ICaLecav
*
* Y(137)   O  Ca channel variable
*
      YDOT(137) =
     *   kco*Y(142) - koc*Y(137)
     * + SFOI1*Kpcb*Y(143) - SFOI1*gammae*Y(137)
     * + SFOI2*alpha*Y(144) - SFOI2*Kpcf*Y(137)
     * - kICaLPKA*Y(89)*Y(137)/(KICaLmPKA+ICaLecav*Y(137))
     * + (alpha/alphap)*kICaLPP*PP1ecav*Y(146)
     * /(KICaLmPP+ICaLecav*Y(146))
*
* Y(138)   C1  Ca channel variable
*
*      YDOT(138) = beta*Y(139)-4.0*alpha*Y(138)
*
* Y(139)   C2  Ca channel variable
*
      YDOT(139) =
     *    4.0*alpha*Y(138) - beta*Y(139)
     *  + 2.0*beta*Y(140) - 3.0*alpha*Y(139)
     * - kICaLPKA*Y(89)*Y(139)/(KICaLmPKA+ICaLecav*Y(139))
     * + (alphap*alphap*kcop/(alpha*alpha*kco))*kICaLPP*PP1ecav*Y(148)
     * /(KICaLmPP+ICaLecav*Y(148))
*
* Y(140)   C3  Ca channel variable
*
      YDOT(140) =
     *    3.0*alpha*Y(139) - 2.0*beta*Y(140)
     *  + 3.0*beta*Y(141) - 2.0*alpha*Y(140)
     * - kICaLPKA*Y(89)*Y(140)/(KICaLmPKA+ICaLecav*Y(140))
     * + (alphap*kcop/(alpha*kco))*kICaLPP*PP1ecav*Y(149)
     * /(KICaLmPP+ICaLecav*Y(149))
*
* Y(141)   C4  Ca channel variable
*
      YDOT(141) =
     *    2.0*alpha*Y(140) - 3.0*beta*Y(141)
     *  + 4.0*beta*Y(142) - alpha*Y(141)
     *  + 4.0*SFC4I1*Kpcb*beta*Y(143) 
     *  - SFC4I1*alpha*gammae*(kco/koc)*Y(141)
     *  + 4.0*SFC4I2*beta*Y(144) - SFC4I2*Kpcf*(kco/koc)*Y(141)
     *  + 4.0*SFC4I3*beta*Kpcb*Y(145) 
     *  - SFC4I3*gammae*Kpcf*(kco/koc)*Y(141)
     * - kICaLPKA*Y(89)*Y(141)/(KICaLmPKA+ICaLecav*Y(141))
     * + (kcop/kco)*kICaLPP*PP1ecav*Y(150)
     * /(KICaLmPP+ICaLecav*Y(150))
*
* Y(142) CP Ca channel variable
*
      YDOT(142) =
     *    alpha*Y(141) - 4.0*beta*Y(142)
     *  + koc*Y(137) - kco*Y(142)
     * - kICaLPKA*Y(89)*Y(142)/(KICaLmPKA+ICaLecav*Y(142))
     * + (alpha*kcop/(alphap*kco))*kICaLPP*PP1ecav*Y(151)
     * /(KICaLmPP+ICaLecav*Y(151))
*
* Y(143)   I1  Ca channel variable
*
      YDOT(143) =
     *    SFOI1*gammae*Y(137) - SFOI1*Kpcb*Y(143)
     *  + SFI1I3*alpha*Y(145) - SFI1I3*Kpcf*Y(143)
     *  + SFC4I1*alpha*gammae*(kco/koc)*Y(141) 
     *  - 4.0*SFC4I1*beta*Kpcb*Y(143)
     * - kICaLPKA*Y(89)*Y(143)/(KICaLmPKA+ICaLecav*Y(143))
     * + (alpha/alphap)*kICaLPP*PP1ecav*Y(152)
     * /(KICaLmPP+ICaLecav*Y(152))
*
* Y(144)   I2  Ca channel variable
*
      YDOT(144) =
     *    SFOI2*Kpcf*Y(137) - SFOI2*alpha*Y(144)
     *  + SFI2I3*Kpcb*Y(145) - SFI2I3*gammae*Y(144)
     *  + SFC4I2*Kpcf*(kco/koc)*Y(141) - 4.0*SFC4I2*beta*Y(144)
     * - kICaLPKA*Y(89)*Y(144)/(KICaLmPKA+ICaLecav*Y(144))
     * + kICaLPP*PP1ecav*Y(153)/(KICaLmPP+ICaLecav*Y(153))
*
* Y(145)   I3  Ca channel variable
*
      YDOT(145) =
     *    SFI1I3*Kpcf*Y(143) - SFI1I3*alpha*Y(145)
     *  + SFI2I3*gammae*Y(144) - SFI2I3*Kpcb*Y(145)
     *  + SFC4I3*gammae*Kpcf*(kco/koc)*Y(141) 
     *  - 4.0*SFC4I3*beta*Kpcb*Y(145)
     * - kICaLPKA*Y(89)*Y(145)/(KICaLmPKA+ICaLecav*Y(145))
     * + kICaLPP*PP1ecav*Y(154)/(KICaLmPP+ICaLecav*Y(154))
*
* Y(146) OP Ca channel variable
*
      YDOT(146) =
     *   kcop*Y(151) - koc*Y(146)
     * + SFOI1*Kpcb*Y(152) - SFOI1*gammae*Y(146)
     * + SFOI2*alphap*Y(153) - SFOI2*Kpcf*Y(146)
     * + kICaLPKA*Y(89)*Y(137)/(KICaLmPKA+ICaLecav*Y(137))  
     * - (alpha/alphap)*kICaLPP*PP1ecav*Y(146)
     * /(KICaLmPP+ICaLecav*Y(146))
*
* Y(147) C1P Ca channel variable
*
      YDOT(147) = beta*Y(148)-4.0*alphap*Y(147)
     * + kICaLPKA*Y(89)*Y(138)/(KICaLmPKA+ICaLecav*Y(138))
     * - (alphap*alphap*alphap*kcop/(alpha*alpha*alpha*kco))
     * *kICaLPP*PP1ecav*Y(147)/(KICaLmPP+ICaLecav*Y(147))
*
* Y(148) C2P Ca channel variable
*
      YDOT(148) =
     *    4.0*alphap*Y(147) - beta*Y(148)
     *  + 2.0*beta*Y(149) - 3.0*alphap*Y(148)
     * + kICaLPKA*Y(89)*Y(139)/(KICaLmPKA+ICaLecav*Y(139))
     * - (alphap*alphap*kcop/(alpha*alpha*kco))*kICaLPP*PP1ecav*Y(148)
     * /(KICaLmPP+ICaLecav*Y(148))
*
* Y(149) C3P Ca channel variable
*
      YDOT(149) =
     *    3.0*alphap*Y(148) - 2.0*beta*Y(149)
     *  + 3.0*beta*Y(150) - 2.0*alphap*Y(149)
     * + kICaLPKA*Y(89)*Y(140)/(KICaLmPKA+ICaLecav*Y(140))
     * - (alphap*kcop/(alpha*kco))*kICaLPP*PP1ecav*Y(149)
     * /(KICaLmPP+ICaLecav*Y(149))
*
* Y(150) C4P Ca channel variable
*
      YDOT(150) =
     *    2.0*alphap*Y(149) - 3.0*beta*Y(150)
     *  + 4.0*beta*Y(151) - alphap*Y(150)
     *  + 4.0*SFC4I1*Kpcb*beta*Y(152)
     *  - SFC4I1*alphap*gammae*(kcop/koc)*Y(150)
     *  + 4.0*SFC4I2*beta*Y(153) - SFC4I2*Kpcf*(kcop/koc)*Y(150)
     *  + 4.0*SFC4I3*beta*Kpcb*Y(154)
     *  - SFC4I3*gammae*Kpcf*(kcop/koc)*Y(150)
     * + kICaLPKA*Y(89)*Y(141)/(KICaLmPKA+ICaLecav*Y(141))
     * - (kcop/kco)*kICaLPP*PP1ecav*Y(150)
     * /(KICaLmPP+ICaLecav*Y(150))
*
* Y(151) CPP Ca channel variable
*
      YDOT(151) =
     *    alphap*Y(150) - 4.0*beta*Y(151)
     *  + koc*Y(146) - kcop*Y(151)
     * + kICaLPKA*Y(89)*Y(142)/(KICaLmPKA+ICaLecav*Y(142))
     * - (alpha*kcop/(alphap*kco))*kICaLPP*PP1ecav*Y(151)
     * /(KICaLmPP+ICaLecav*Y(151))
*
* Y(152) I1P Ca channel variable
*
      YDOT(152) =
     *    SFOI1*gammae*Y(146) - SFOI1*Kpcb*Y(152)
     *  + SFI1I3*alphap*Y(154) - SFI1I3*Kpcf*Y(152)
     *  + SFC4I1*alphap*gammae*(kcop/koc)*Y(150)
     *  - 4.0*SFC4I1*beta*Kpcb*Y(152)
     * + kICaLPKA*Y(89)*Y(143)/(KICaLmPKA+ICaLecav*Y(143))
     * - (alpha/alphap)*kICaLPP*PP1ecav*Y(152)
     * /(KICaLmPP+ICaLecav*Y(152))
*
* Y(153) I2P Ca channel variable
*
      YDOT(153) =
     *    SFOI2*Kpcf*Y(146) - SFOI2*alphap*Y(153)
     *  + SFI2I3*Kpcb*Y(154) - SFI2I3*gammae*Y(153)
     *  + SFC4I2*Kpcf*(kcop/koc)*Y(150) - 4.0*SFC4I2*beta*Y(153)
     * + kICaLPKA*Y(89)*Y(144)/(KICaLmPKA+ICaLecav*Y(144))
     * - kICaLPP*PP1ecav*Y(153)/(KICaLmPP+ICaLecav*Y(153))
*
* Y(154) I3P Ca channel variable
*
      YDOT(154) =
     *    SFI1I3*Kpcf*Y(152) - SFI1I3*alphap*Y(154)
     *  + SFI2I3*gammae*Y(153) - SFI2I3*Kpcb*Y(154)
     *  + SFC4I3*gammae*Kpcf*(kcop/koc)*Y(150)
     *  - 4.0*SFC4I3*beta*Kpcb*Y(154)
     * + kICaLPKA*Y(89)*Y(145)/(KICaLmPKA+ICaLecav*Y(145))
     * - kICaLPP*PP1ecav*Y(154)/(KICaLmPP+ICaLecav*Y(154))
*
* Y(16)-Y(19)  RyR channel states
*
      temp5 = Y(7)*Y(7)*Y(7)
      temp6 = Y(7)*Y(7)*Y(7)*Y(7)
*
* Y(16) RyR C1 channel variable
*
      YDOT(16)=-Kap*temp6*Y(16)+Kam*Y(18)
     * - kRyRPKA*Y(89)*Y(16)/(KRyRmPKA+RyRecav*Y(16))
     * + kRyRPP*PP1ecav*Y(116)/(KRyRmPP+RyRecav*Y(116))
*
* Y(17) RyR C2 channel variable
*
      YDOT(17)= Kcp*Y(18)-Kcm*Y(17)
     * - kRyRPKA*Y(89)*Y(17)*f_RyR/(KRyRmPKA+RyRecav*Y(17))
     * + ((Kap*Kamp*Kcp*Kcmp)/(Kapp*Kam*Kcpp*Kcm))*f_RyR
     * * kRyRPP*PP1ecav*Y(117)/(KRyRmPP+RyRecav*Y(117))
*
* Y(18) RyR O1 channel variable
*
      YDOT(18)= Kap*temp6*Y(16)-Kam*Y(18)
     *         -Kbp*temp5*Y(18)+Kbm*Y(19)
     *         -Kcp*Y(18)+Kcm*Y(17)
     * - kRyRPKA*Y(89)*Y(18)*f_RyR/(KRyRmPKA+RyRecav*Y(18))
     * + ((Kap*Kamp)/(Kapp*Kam))*f_RyR*kRyRPP*PP1ecav*Y(118)
     * /(KRyRmPP+RyRecav*Y(118))
*
* Y(19) RyR O2 channel variable
*
      YDOT(19)= Kbp*temp5*Y(18)-Kbm*Y(19)
     * - kRyRPKA*Y(89)*Y(19)*f_RyR/(KRyRmPKA+RyRecav*Y(19))
     * + ((Kap*Kamp*Kbp*Kbmp)/(Kapp*Kam*Kbpp*Kbm))*f_RyR
     * * kRyRPP*PP1ecav*Y(119)/(KRyRmPP+RyRecav*Y(119))
*
* Y(116) RyR C1P channel variable
*
      YDOT(116)=-Kapp*temp6*Y(116)+Kamp*Y(118)
     * + kRyRPKA*Y(89)*Y(16)/(KRyRmPKA+RyRecav*Y(16))
     * - kRyRPP*PP1ecav*Y(116)/(KRyRmPP+RyRecav*Y(116))
*
* Y(117) RyR C2P channel variable
*
      YDOT(117)= Kcpp*Y(118)-Kcmp*Y(117)
     * + kRyRPKA*Y(89)*Y(17)*f_RyR/(KRyRmPKA+RyRecav*Y(17))
     * - ((Kap*Kamp*Kcp*Kcmp)/(Kapp*Kam*Kcpp*Kcm))*f_RyR
     * * kRyRPP*PP1ecav*Y(117)/(KRyRmPP+RyRecav*Y(117))
*
* Y(118) RyR O1P channel variable
*
      YDOT(118)= Kapp*temp6*Y(116)-Kamp*Y(118)
     *          -Kbpp*temp5*Y(118)+Kbmp*Y(119)
     *          -Kcpp*Y(118)+Kcmp*Y(117)
     * + kRyRPKA*Y(89)*Y(18)*f_RyR/(KRyRmPKA+RyRecav*Y(18))
     * - ((Kap*Kamp)/(Kapp*Kam))*f_RyR*kRyRPP*PP1ecav*Y(118)
     * /(KRyRmPP+RyRecav*Y(118))
*
* Y(119) RyR O2P channel variable
*
      YDOT(119)= Kbpp*temp5*Y(118)-Kbmp*Y(119)
     * + kRyRPKA*Y(89)*Y(19)*f_RyR/(KRyRmPKA+RyRecav*Y(19))
     * - ((Kap*Kamp*Kbp*Kbmp)/(Kapp*Kam*Kbpp*Kbm))*f_RyR
     * * kRyRPP*PP1ecav*Y(119)/(KRyRmPP+RyRecav*Y(119))
*
*  Y(20)-Y(22),Y(37)-Y(42) Na fast current (Clancy-Rudy, Circulation, 2002)
*
      Va=Y(1)-2.5
      Vi=Y(1)+7.0
      alp11 = 3.802/(0.1027*exp(-Va/17.0)+0.20*exp(-Va/150.))
      alp12 = 3.802/(0.1027*exp(-Va/15.0)+0.23*exp(-Va/150.))
      alp13 = 3.802/(0.1027*exp(-Va/12.0)+0.25*exp(-Va/150.))
      bet11 = 0.1917*exp(-Va/20.3)
      bet12 = 0.20*exp(-(Va-5.0)/20.3)
      bet13 = 0.22*exp(-(Va-10.0)/20.3)
      alp3  = 7.0e-7*exp(-Vi/7.7)
      bet3  = (0.0084+0.00002*Vi)
      alp2  = 1.0/(0.188495*exp(-Vi/16.6)+0.393956)
      bet2  = alp13*alp2*alp3/(bet13*bet3)
      alp4  = alp2/100.
      bet4  = alp3
      alp5  = alp2/9.5e4
      bet5  = alp3/50.0
      Vap=Y(1)-2.5
      Vip=Y(1)+7.0
      alp11p = 3.802/(0.1027*exp(-Vap/17.0)+0.20*exp(-Vap/150.))
      alp12p = 3.802/(0.1027*exp(-Vap/15.0)+0.23*exp(-Vap/150.))
      alp13p = 3.802/(0.1027*exp(-Vap/12.0)+0.25*exp(-Vap/150.))
      bet2p  = alp13p*alp2*alp3/(bet13*bet3)
*
* Y(20) INa CNa3 channel variable
*
*      YDOT(20) = bet11*Y(21) - alp11*Y(20)
*     *         + alp3*Y(42)  - bet3*Y(20)
*     * - kINAPKA*Y(83)*Y(20)/(KINamPKA+Y(20))
*     * + (alp11p*alp12p/(alp11*alp12))*kINaPP*PPcav*Y(120)
*     * /(KINamPP+Y(120))
*
* Y(21) INa CNa2 channel variable
*
      YDOT(21) = alp11*Y(20) - bet11*Y(21)
     *         + bet12*Y(22) - alp12*Y(21)
     *         + alp3*Y(41)  - bet3*Y(21)
     * - kINAPKA*Y(83)*Y(21)/(KINamPKA+Y(21))
     * + (alp12p/alp12)*kINaPP*PPcav*Y(121)/(KINamPP+Y(121))
*
* Y(22) INa CNa1 channel variable
*
      YDOT(22) = alp12*Y(21) - bet12*Y(22)
     *         + bet13*Y(37) - alp13*Y(22)
     *         + alp3*Y(38)  - bet3*Y(22)
     * - kINAPKA*Y(83)*Y(22)/(KINamPKA+Y(22))
     * + kINaPP*PPcav*Y(122)/(KINamPP+Y(122))
*
* Y(37) INa ONa channel variable
*
      YDOT(37) = alp13*Y(22) - bet13*Y(37)
     *         + bet2*Y(38)  - alp2*Y(37)
     * - kINAPKA*Y(83)*Y(37)/(KINamPKA+Y(37))
     * + (alp13/alp13p)*kINaPP*PPcav*Y(123)/(KINamPP+Y(123))
*
* Y(38) INa IFNa channel variable
*
      YDOT(38) = alp2*Y(37)  - bet2*Y(38)
     *         + bet3*Y(22)  - alp3*Y(38)
     *         + bet4*Y(39)  - alp4*Y(38)
     *         + alp12*Y(41) - bet12*Y(38)
     * - kINAPKA*Y(83)*Y(38)/(KINamPKA+Y(38))
     * + kINaPP*PPcav*Y(124)/(KINamPP+Y(124))
*
* Y(39) INa I1Na channel variable
*
      YDOT(39) = alp4*Y(38)  - bet4*Y(39)
     *         + bet5*Y(40)  - alp5*Y(39)
     * - kINAPKA*Y(83)*Y(39)/(KINamPKA+Y(39))
     * + kINaPP*PPcav*Y(125)/(KINamPP+Y(125))
*
* Y(40) INa I2Na channel variable
*
      YDOT(40) = alp5*Y(39)  - bet5*Y(40)
     * - kINAPKA*Y(83)*Y(40)/(KINamPKA+Y(40))
     * + kINaPP*PPcav*Y(126)/(KINamPP+Y(126))
*
* Y(41) INa ICNa2 channel variable
*
      YDOT(41) = alp11*Y(42) - bet11*Y(41)
     *         + bet12*Y(38) - alp12*Y(41)
     *         + bet3*Y(21)  - alp3*Y(41)
     * - kINAPKA*Y(83)*Y(41)/(KINamPKA+Y(41))
     * + (alp12p/alp12)*kINaPP*PPcav*Y(127)/(KINamPP+Y(127))
*
* Y(42) INa ICNa3 channel variable
*
      YDOT(42) = bet11*Y(41) - alp11*Y(42)
     *         + bet3*Y(20)  - alp3*Y(42)
     * - kINAPKA*Y(83)*Y(42)/(KINamPKA+Y(42))
     * + (alp11p*alp12p/(alp11*alp12))*kINaPP*PPcav*Y(128)
     * /(KINamPP+Y(128))
*
* Y(120) INa CNa3p channel variable
*
      YDOT(120) = bet11*Y(121) - alp11p*Y(120)
     *          + alp3*Y(128)  - bet3*Y(120)
     * + kINAPKA*Y(83)*Y(20)/(KINamPKA+Y(20))
     * - (alp11p*alp12p/(alp11*alp12))*kINaPP*PPcav*Y(120)
     * /(KINamPP+Y(120))
*
* Y(121) INa CNa2p channel variable
*
      YDOT(121) = alp11p*Y(120) - bet11*Y(121)
     *          + bet12*Y(122)  - alp12p*Y(121)
     *          + alp3*Y(127)   - bet3*Y(121) 
     * + kINAPKA*Y(83)*Y(21)/(KINamPKA+Y(21))
     * - (alp12p/alp12)*kINaPP*PPcav*Y(121)/(KINamPP+Y(121))
*
* Y(122) INa CNa1p channel variable
*
      YDOT(122) = alp12p*Y(121) - bet12*Y(122)
     *          + bet13*Y(123)  - alp13p*Y(122)
     *          + alp3*Y(124)   - bet3*Y(122)
     * + kINAPKA*Y(83)*Y(22)/(KINamPKA+Y(22))
     * - kINaPP*PPcav*Y(122)/(KINamPP+Y(122))
*
* Y(123) INa ONap channel variable
*
      YDOT(123) = alp13p*Y(122) - bet13*Y(123)
     *          + bet2p*Y(124)  - alp2*Y(123)
     * + kINAPKA*Y(83)*Y(37)/(KINamPKA+Y(37))
     * - (alp13/alp13p)*kINaPP*PPcav*Y(123)/(KINamPP+Y(123))
*
* Y(124) INa IFNap channel variable
*
      YDOT(124) = alp2*Y(123)   - bet2p*Y(124)
     *          + bet3*Y(122)   - alp3*Y(124)
     *          + bet4*Y(125)   - alp4*Y(124)
     *          + alp12p*Y(127) - bet12*Y(124)
     * + kINAPKA*Y(83)*Y(38)/(KINamPKA+Y(38))
     * - kINaPP*PPcav*Y(124)/(KINamPP+Y(124))
*
* Y(125) INa I1Nap channel variable
*
      YDOT(125) = alp4*Y(124)  - bet4*Y(125)
     *          + bet5*Y(126)  - alp5*Y(125)
     * + kINAPKA*Y(83)*Y(39)/(KINamPKA+Y(39))
     * - kINaPP*PPcav*Y(125)/(KINamPP+Y(125))
*
* Y(126) INa I2Nap channel variable
*
      YDOT(126) = alp5*Y(125)  - bet5*Y(126)
     * + kINAPKA*Y(83)*Y(40)/(KINamPKA+Y(40))
     * - kINaPP*PPcav*Y(126)/(KINamPP+Y(126))
*
* Y(127) INa ICNa2p channel variable
*
      YDOT(127) = alp11p*Y(128) - bet11*Y(127)
     *          + bet12*Y(124)  - alp12p*Y(127)
     *          + bet3*Y(121)   - alp3*Y(127)
     * + kINAPKA*Y(83)*Y(41)/(KINamPKA+Y(41))
     * - (alp12p/alp12)*kINaPP*PPcav*Y(127)/(KINamPP+Y(127))
*
* Y(128) INa ICNa3p channel variable
*
      YDOT(128) = bet11*Y(127) - alp11p*Y(128)
     *          + bet3*Y(120)  - alp3*Y(128)
     * + kINAPKA*Y(83)*Y(42)/(KINamPKA+Y(42))
     * - (alp11p*alp12p/(alp11*alp12))*kINaPP*PPcav*Y(128)
     * /(KINamPP+Y(128))
*
*  Y(23)  Na intracellular concentration
*
      YDOT(23) = -1.0*(Ina+Inab+3.0*(Inaca+Inak))*Acap/(Vmyo*F)
*
*  Y(24)  K  intracellular concentration
*
      YDOT(24) = -1.0*(Ikto+Ik1+Iks+Ikur+Ikr+ICaK
     *  -2.0*Inak-Istim)*Acap/(Vmyo*F)
*
*  Y(25),Y(26)  ato and ito gating variables for Ikto
*
      alp25 = 0.04516*exp(0.03577*(Y(1)+33.0))*4.0
      bet25 = 0.0989*exp(-0.06237*(Y(1)+33.0))*4.0
       temp7 = 0.0019*exp((Y(1)+15.5)/(-1.0*dito))
       temp8 = 0.067083*exp((Y(1)+15.5+20.0)/(-1.0*dito))
      alp26 = 0.08*temp7/(1.0+temp8)
       temp9  = 0.0019*exp((Y(1)+15.5+20.0)/dito)
       temp10 = 0.051335*exp((Y(1)+15.5+20.0)/dito)
      bet26 = 0.5*temp9/(1.0+temp10)
*
      YDOT(25) = alp25*(1.0-Y(25)) - bet25*Y(25)
      YDOT(26) = alp26*(1.0-Y(26)) - bet26*Y(26)
*
*  Y(135),Y(136) ato and ito gating variables for Ikto,phosphorylated
*
      alp25p = 0.04516*exp(0.03577*(Y(1)+30.0-13.0))*4.0
      bet25p = 0.0989*exp(-0.06237*(Y(1)+30.0-13.0))*4.0
       temp7p = 0.0019*exp((Y(1)+13.5-6.0)/(-1.0*dito))
       temp8p = 0.067083*exp((Y(1)+13.5+20.0-6.0)/(-1.0*dito))
      alp26p = 0.08*temp7p/(1.0+temp8p)
       temp9p  = 0.0019*exp((Y(1)+13.5+20.0-6.0)/dito)
       temp10p = 0.051335*exp((Y(1)+13.5+20.0-6.0)/dito)
      bet26p = 0.5*temp9p/(1.0+temp10p)
*
      YDOT(135) = alp25p*(1.0-Y(135)) - bet25p*Y(135)
      YDOT(136) = alp26p*(1.0-Y(136)) - bet26p*Y(136)
*
*  Y(27)  nks gating variable for Iks
*
      temp11 = 0.00001444*(Y(1)+26.5)
      alp27  = temp11/(1.0-exp(-0.128*(Y(1)+26.5)))/3.
      bet27  = 0.000286*exp(-0.038*(Y(1)+26.5))/3.
*
      YDOT(27) = alp27*(1.0-Y(27)) - bet27*Y(27)
*
*  Y(28), Y(29)  aur and iur gating variables for Ikur1
*
      ass = 1.0/(1.0+exp((-22.5-Y(1))/7.7))
      iss1 = 1.0/(1.0+exp((45.2+Y(1))/5.7))
*      taua1 = (0.493*exp(-0.0629*Y(1))+2.058)
      taua1 = 6.1/(exp(0.0629*(Y(1)+40.0))
     *      + exp(-0.0629*(Y(1)+40.0)))+2.058
      taui1 = 270.+1050/(1.0+exp((45.2+Y(1))/5.7))
*
      YDOT(28) = (ass-Y(28))/taua1
      YDOT(29) = (iss1-Y(29))/taui1
*
*  Y(35), Y(36)  aur and iur gating variables for Ikur2
*
      iss2 = 1.0/(1.0+exp((45.2+Y(1))/5.7))
*      taua2 = (0.493*exp(-0.0629*Y(1))+2.058)
      taua2 = 6.1/(exp(0.0629*(Y(1)+40.0))
     *      + exp(-0.0629*(Y(1)+40.0)))+2.058
      taui2 = 803.-18./(1.0+exp((45.2+Y(1))/5.7))
*
      YDOT(35) = (ass-Y(35))/taua2
      YDOT(36) = (iss2-Y(36))/taui2
      YDOT(130) = (ass-Y(130))/taua2
      YDOT(131) = (iss2-Y(131))/taui2
      YDOT(132) = kIKurPP*PP1ecav*(1-Y(132))/(KIKurmPP+(1-Y(132)))
     *          - kIKurPKA*Y(89)*Y(132)/(KIKurmPKA+Y(132))
*
*  Y(43) aur gating variable for Ikur3
*
*      taua3 = (39.3*exp(-0.0862*Y(1))+13.17)
      taua3 = 1235.5/(exp(0.0862*(Y(1)+40.0))
     *      + exp(-0.0862*(Y(1)+40.0)))+13.17
      YDOT(43) = (ass-Y(43))/taua3
*
*  Y(184), Y(185) aur and iur gating variables for Ikur4
*
      taui4 = 5334.-4912./(1.0+exp((45.2+Y(1))/5.7))
      YDOT(184) = (ass-Y(184))/taua2
      YDOT(185) = (iss2-Y(185))/taui4
*
* Y(30)-Y(34) HERG channel state variables
*
      ala0 = 0.022348*exp(0.01176*Y(1))
      bea0 = 0.047002*exp(-0.0631*Y(1))
      ala1 = 0.013733*exp(0.038198*Y(1))
      bea1 = 0.0000689*exp(-0.04178*Y(1))
      ali  = 0.090821*exp(0.023391*(Y(1)+5.0))
      bei  = 0.006497*exp(-0.03268*(Y(1)+5.0))
*
      YDOT(30) = bea0*Y(31)-ala0*Y(30)
      YDOT(31) = ala0*Y(30)-bea0*Y(31)+kb*Y(32)-kf*Y(31)
      YDOT(32) = kf*Y(31)-kb*Y(32)+bea1*Y(33)-ala1*Y(32)
      YDOT(33) = ala1*Y(32)-bea1*Y(33)+bei*Y(34)-ali*Y(33)
      YDOT(34) = ali*Y(33)-bei*Y(34)
*
*   Y(44) Pryr Ryanodine receptor modulation factor
*
      YDOT(44) = -t1*Y(44)+t2*Icas*exp((Y(1)+5.0)*(Y(1)+5.0)/(-648.))
*
*   Caveolae domain
*
*   Y(45) beta1 tot concentration phosphorylated by PKA caveolae
*
      YDOT(45) = kPKAp*Y(83)*RB1cavnptot-kPKAm*Y(45)
*   
*   Y(46) beta1 tot concentration phosphorylated by BARK caveolae
*
      YDOT(46) = kBARKp*(LRB1cavnp+LRB1Gscavnp)-kBARKm*Y(46)
*
*   Y(155) beta2 tot concentration phosphorylated by PKA caveolae
*
      YDOT(155) = kPKAp*Y(83)*RB2cavnptot-kPKAm*Y(155)
*
*   Y(156) beta2 tot concentration phosphorylated by BARK caveolae
*
      YDOT(156) = kBARKp*(LRB2cavnp+LRB2Gscavnp)-kBARKm*Y(156)
*
*   Y(47) Gs-alpha with GTP caveolae
*
      YDOT(47) = kact2Gs*RB1Gscavnp+factGsGi*kact2Gs*RB2Gscavnp
     * +kact1Gs*LRB1Gscavnp+factGsGi*kact1Gs*LRB2Gscavnp-khydGs*Y(47)
*
*   Y(48)  G Beta-Gamma caveolae
*
      YDOT(48) = kact2Gs*RB1Gscavnp+factGsGi*kact2Gs*RB2Gscavnp
     *     +kact1Gs*LRB1Gscavnp+factGsGi*kact1Gs*LRB2Gscavnp
     *     +kact2Gi*RB2GicavPKA+kact1Gi*LRB2GicavPKA
     *     -kreasGs*Y(48)*Y(49)-kreasGi*Y(48)*Y(158)
*
*   Y(49)  Gs-alpha with GDP caveolae
*
      YDOT(49) = khydGs*Y(47)-kreasGs*Y(48)*Y(49)
*
*   Y(157) Gi-alpha with GTP caveolae
*
      YDOT(157) = kact2Gi*RB2GicavPKA+kact1Gi*LRB2GicavPKA
     *   -khydGi*Y(157)
*
*   Y(158)  Gi-alpha with GDP caveolae
*
      YDOT(158) = khydGi*Y(157)-kreasGi*Y(48)*Y(158)
*
*   Y(50) beta1 tot concentration phosphorylated by PKA extracaveolae
*
      YDOT(50) = kPKAp*Y(89)*RB1ecavnptot-kPKAm*Y(50)
*
*   Y(51) beta1 tot concentration phosphorylated by BARK extracaveolae
*
      YDOT(51) = kBARKp*(LRB1ecavnp+LRB1Gsecavnp)-kBARKm*Y(51)
*
*   Y(159) beta2 tot concentration phosphorylated by PKA extracaveolae
*
      YDOT(159) = kPKAp*Y(89)*RB2ecavnptot-kPKAm*Y(159)
*
*   Y(160) beta2 tot concentration phosphorylated by BARK extracaveolae
*
      YDOT(160) = kBARKp*(LRB2ecavnp+LRB2Gsecavnp)-kBARKm*Y(160)
*
*   Y(52) Gs-alpha with GTP extracaveolae
*
      YDOT(52) = kact2Gs*RB1Gsecavnp+factGsGi*kact2Gs*RB2Gsecavnp
     *         +kact1Gs*LRB1Gsecavnp+factGsGi*kact1Gs*LRB2Gsecavnp
     *         -khydGs*Y(52)
*
*   Y(53)  G Beta-Gamma extracaveolae
*
      YDOT(53) = kact2Gs*RB1Gsecavnp+factGsGi*kact2Gs*RB2Gsecavnp
     *     +kact1Gs*LRB1Gsecavnp+factGsGi*kact1Gs*LRB2Gsecavnp
     *     +kact2Gi*RB2GiecavPKA+kact1Gi*LRB2GiecavPKA
     *     -kreasGs*Y(53)*Y(54)-kreasGi*Y(53)*Y(162)
*
*   Y(54)  Gs-alpha with GDP extracaveolae
*
      YDOT(54) = khydGs*Y(52)-kreasGs*Y(53)*Y(54)
*
*   Y(161) Gi-alpha with GTP extracaveolae
*
      YDOT(161) = kact2Gi*RB2GiecavPKA+kact1Gi*LRB2GiecavPKA
     *   -khydGi*Y(161)
*
*   Y(162)  Gi-alpha with GDP extracaveolae
*
      YDOT(162) = khydGi*Y(161)-kreasGi*Y(53)*Y(162)
*
*   Y(55) beta1 tot concentration phosphorylated by PKA cytosol
*
      YDOT(55) = kPKAp*Y(95)*RB1cytnptot-kPKAm*Y(55)
*
*   Y(56) beta1 tot concentration phosphorylated by BARK
*
      YDOT(56) = kBARKp*(LRB1cytnp+LRB1Gscytnp)-kBARKm*Y(56)
*
*   Y(57) Gs-alpha with GTP
*
      YDOT(57) = kact2Gs*RB1Gscytnp+kact1Gs*LRB1Gscytnp
     *         -khydGs*Y(57)
*
*   Y(58)  Gs Beta-Gamma
*
      YDOT(58) = kact2Gs*RB1Gscytnp+kact1Gs*LRB1Gscytnp
     *         - kreasGs*Y(58)*Y(59)
*
*   Y(59)  Gs-alpha with GDP
*
      YDOT(59) = khydGs*Y(57)-kreasGs*Y(58)*Y(59)
*
*   Y(60) cAMP from AC56 in ceveolae
*
      YDOT(60) = kcavAC56*AC56cav*ATP/(Kmatp+ATP)
*
*   Y(61) cAMP from AC47 in extracaveolae
*
      YDOT(61) = kecavAC47*AC47ecav*ATP/(Kmatp+ATP)
*
*   Y(62) cAMP from AC56 in cytosol
*
      YDOT(62) = kcytAC56*AC56cyt*ATP/(Kmatp+ATP)
*
*   Y(63) cAMP from AC47 in cytosol
*
      YDOT(63) = kcytAC47*AC47cyt*ATP/(Kmatp+ATP)
*
*   Y(64) PDE3 caveolar phosphorylated
*
      YDOT(64) = kfpdep*Y(83)*(PDE3cavtot-Y(64))-kbpdep*Y(64)
*
*   Y(65) PDE4 caveolar phosphorylated
*
      YDOT(65) = kfpdep*Y(83)*(PDE4cavtot-Y(65))-kbpdep*Y(65)
*
*   Y(66) cAMP change due to PDE2 caveolar domain
*
      YDOT(66) = (kpde2*PDE2cavtot*Y(97))/(Kmpde2+Y(97))
*
*   Y(67) cAMP change due to PDE3 caveolar domain
*
      YDOT(67) = (kpde3*(PDE3cavtot-Y(64))*Y(97)+deltakpde34*kpde3
     *           *Y(64)*Y(97))/(Kmpde3+Y(97))
*
*   Y(68) cAMP change due to PDE4 caveolar domain
*
      YDOT(68) = (kpde4*(PDE4cavtot-Y(65))*Y(97)+deltakpde34*kpde4
     *           *Y(65)*Y(97))/(Kmpde4+Y(97))
*
*   Y(69) PDE3 extracaveolar phosphorylated
*
*      YDOT(69) = kfpdep*Y(89)*(PDE3ecavtot-Y(69))-kbpdep*Y(69)
       YDOT(69) = 0.0
*
*   Y(70) PDE4 extracaveolar phosphorylated
*
      YDOT(70) = kfpdep*Y(89)*(PDE4ecavtot-Y(70))-kbpdep*Y(70)
*
*   Y(71) cAMP change due to PDE2 extracaveolar domain
*
      YDOT(71) = (kpde2*PDE2ecavtot*Y(98))/(Kmpde2+Y(98))
*
*   Y(72) cAMP change due to PDE3 extracaveolar domain
*
      YDOT(72) = (kpde3*(PDE3ecavtot-Y(69))*Y(98)+deltakpde34*kpde3
     *           *Y(69)*Y(98))/(Kmpde3+Y(98))
*
*   Y(73) cAMP change due to PDE4 extracaveolar domain
*
      YDOT(73) = (kpde4*(PDE4ecavtot-Y(70))*Y(98)+deltakpde34*kpde4
     *           *Y(70)*Y(98))/(Kmpde4+Y(98))
*
*   Y(74) PDE3 cytosol phosphorylated
*
      YDOT(74) = kfpdep*Y(95)*(PDE3cyttot-Y(74))-kbpdep*Y(74)
*
*   Y(75) PDE4 cytosol phosphorylated
*
      YDOT(75) = kfpdep*Y(95)*(PDE4cyttot-Y(75))-kbpdep*Y(75)
*
*   Y(76) cAMP change due to PDE2 cytosol domain
*
      YDOT(76) = (kpde2*PDE2cyttot*Y(99))/(Kmpde2+Y(99))
*
*   Y(77) cAMP change due to PDE3 cytosol domain
*
      YDOT(77) = (kpde3*(PDE3cyttot-Y(74))*Y(99)+deltakpde34*kpde3
     *           *Y(74)*Y(99))/(Kmpde3+Y(99))
*
*   Y(78) cAMP change due to PDE4 cytosol domain
*
      YDOT(78) = (kpde4*(PDE4cyttot-Y(75))*Y(99)+deltakpde34*kpde4
     *           *Y(75)*Y(99))/(Kmpde4+Y(99))
*
*   Y(79) cAMP binding to PKA caveolar
*
      YDOT(79) = -kpkaiif1*RCcavf*Y(97)+kpkaiib1*Y(80)
     *         -kpkaiif2*Y(80)*Y(97)+kpkaiib2*Y(81)
*      write(51,44) time,-kpkaiif1*RCcavf*Y(97),kpkaiib1*Y(80),
*     *                  -kpkaiif2*Y(80)*Y(97),kpkaiib2*Y(81)
*
*   Y(80) RC with bound one cAMP in caveolae
*
      YDOT(80) = kpkaiif1*RCcavf*Y(97)-kpkaiib1*Y(80)
     *         -kpkaiif2*Y(80)*Y(97)+kpkaiib2*Y(81)
*
*   Y(81) RC with bound two cAMP in caveolae
*
      YDOT(81) = kpkaiif2*Y(80)*Y(97)-(kpkaiib2+kpkaiif3)*Y(81)
     *         +kpkaiib3*Y(82)*Y(83)
*
*   Y(82) R with bound two cAMP in caveolae
*
      YDOT(82) = kpkaiif3*Y(81)-kpkaiib3*Y(82)*Y(83)
*
*   Y(83) concentration of catalytic subunit of PKA in caveolae
*
      YDOT(83) = kpkaiif3*Y(81)-kpkaiib3*Y(82)*Y(83)+kpkib*Y(84)
     *         -kpkif*PKIcavf*Y(83)
*
*   Y(84) concentration of PKI bound to catalytic subunit of PKA cav
*
      YDOT(84) = -kpkib*Y(84)+kpkif*PKIcavf*Y(83)
*
*   Extracaveolar domain
*
*   Y(85) cAMP binding to PKA extracaveolar
*
      YDOT(85) = -kpkaiif1*RCecavf*Y(98)+kpkaiib1*Y(86)
     *         -kpkaiif2*Y(86)*Y(98)+kpkaiib2*Y(87)
*
*   Y(86) RC with bound one cAMP in extracaveolae
*
      YDOT(86) = kpkaiif1*RCecavf*Y(98)-kpkaiib1*Y(86)
     *         -kpkaiif2*Y(86)*Y(98)+kpkaiib2*Y(87)
*
*   Y(87) RC with bound two cAMP in extracaveolae
*
      YDOT(87) = kpkaiif2*Y(86)*Y(98)-(kpkaiib2+kpkaiif3)*Y(87)
     *         +kpkaiib3*Y(88)*Y(89)
*
*   Y(88) R with bound two cAMP in extracaveolae
*
      YDOT(88) = kpkaiif3*Y(87)-kpkaiib3*Y(88)*Y(89)
*
*   Y(89) concentration of catalytic subunit of PKA in extracaveolae
*
      YDOT(89) = kpkaiif3*Y(87)-kpkaiib3*Y(88)*Y(89)+kpkib*Y(90)
     *         -kpkif*PKIecavf*Y(89)
*
*   Y(90) concentration of PKI bound to catalytic subunit of PKA ecav
*
      YDOT(90) = -kpkib*Y(90)+kpkif*PKIecavf*Y(89)
*
*   Cytosol domain
*
*   Y(91) cAMP binding to PKA cytosol
*
      YDOT(91) = -kpkaif1*RCcytf*Y(99)+kpkaib1*Y(92)
     *         -kpkaif2*Y(92)*Y(99)+kpkaib2*Y(93)
*
*   Y(92) RC with bound one cAMP in cytosol
*
      YDOT(92) = kpkaif1*RCcytf*Y(99)-kpkaib1*Y(92)
     *         -kpkaif2*Y(92)*Y(99)+kpkaib2*Y(93)
*
*   Y(93) RC with bound two cAMP in cytosol
*
      YDOT(93) = kpkaif2*Y(92)*Y(99)-(kpkaib2+kpkaif3)*Y(93)
     *         +kpkaib3*Y(94)*Y(95)
*
*   Y(94) R with bound two cAMP in cytosol
*
      YDOT(94) = kpkaif3*Y(93)-kpkaib3*Y(94)*Y(95)
*
*   Y(95) concentration of catalytic subunit of PKA in cytosol
*
      YDOT(95) = kpkaif3*Y(93)-kpkaib3*Y(94)*Y(95)+kpkib*Y(96)
     *         -kpkif*PKIcytf*Y(95)
*
*   Y(96) concentration of PKI bound to catalytic subunit of PKA cyt
*
      YDOT(96) = -kpkib*Y(96)+kpkif*PKIcytf*Y(95)
*
*   Y(101) cAMP change due to PDE8 caveolar domain
*
      YDOT(101) = (kpde8*PDE8cavtot*Y(97))/(Kmpde8+Y(97))
*
*   Y(97) total change of cAMP in ceveolae
*
      YDOT(97) =
     *      YDOT(79)+YDOT(60)-YDOT(66)-YDOT(67)-YDOT(68)-YDOT(101)
     *         - jcavecav*(Y(97)-Y(98))/Vcav
     *         - jcavcyt*(Y(97)-Y(99))/Vcav
      fluxcavecav=jcavecav*(Y(97)-Y(98))
      fluxcavcyt =jcavcyt*(Y(97)-Y(99))
*
*   Y(102) cAMP change due to PDE8 extracaveolar domain
*
      YDOT(102) = (kpde8*PDE8ecavtot*Y(98))/(Kmpde8+Y(98))
*
*   Y(98) total change of cAMP in extracaveolae
*
      YDOT(98) =
     *      YDOT(85)+YDOT(61)-YDOT(71)-YDOT(72)-YDOT(73)-YDOT(102)
     *         - jcavecav*(Y(98)-Y(97))/Vecav
     *         - jecavcyt*(Y(98)-Y(99))/Vecav
      fluxecavcav=jcavecav*(Y(98)-Y(97))
      fluxecavcyt=jecavcyt*(Y(98)-Y(99))
*
*   Y(103) cAMP change due to PDE8 cytosolic domain
*
      YDOT(103) = (kpde8*PDE8cyttot*Y(99))/(Kmpde8+Y(99))
*
*   Y(99) total change of cAMP in cytosol
*
      YDOT(99) = YDOT(91)+YDOT(62)+YDOT(63)
     *         - YDOT(76)-YDOT(77)-YDOT(78)-YDOT(103)
     *         - jcavcyt*(Y(99)-Y(97))/Vcyt
     *         - jecavcyt*(Y(99)-Y(98))/Vcyt
      fluxcytcav=jcavcyt*(Y(99)-Y(97))
      fluxcytecav=jecavcyt*(Y(99)-Y(98))
*
*    Y(100) Inhibitor 1 cytosol phosphorylated
*
      YDOT(100)  = kpkainhib1*Y(95)*inhib1cytf
     *             /(Kmpkainhib1+inhib1cytf)-kpp2ainhib1pp2acyt
     *             *Y(100)/(Kmpp2ainhib1+Y(100))
*
*    Y(104) PLB fraction phosphorylated
*
      YDOT(104) = kPLBPKA*Y(95)*(1-Y(104))/(KPLBmPKA+(1-Y(104)))
     *          - kPLBPP1*pp1cytf*Y(104)/(KPLBmPP1+Y(104))
*
*    Y(105) TnI fraction phosphorylated
*
      YDOT(105) = kTnIPKA*Y(95)*(1-Y(105))/(KTnImPKA+(1-Y(105)))
     *          - kTnIPP2A*PP2Acyt*Y(105)/(KTnImPP2A+Y(105))
*
*    Y(129) INaK fraction phosphorylated
*
      YDOT(129) = kINaKPKA*Y(83)*(1-Y(129))/(KINaKmPKA+(1-Y(129)))
     *          - kINaKPP*PPcav*Y(129)/(KINaKmPP+Y(129))
*
*    Y(133) IK1 fraction nonphophorylated
*
      YDOT(133) = kIK1PP*PP1ecav*(1-Y(133))/(KIK1mPP+(1-Y(133)))
     *          - kIK1PKA*Y(89)*Y(133)/(KIK1mPKA+Y(133))
*
*    Y(134) IKto fraction phophorylated
*
      YDOT(134) = kIKtoPKA*Y(89)*(1-Y(134))/(KIKtomPKA+(1-Y(134)))
     *          - kIKtoPP*PP1ecav*Y(134)/(KIKtomPP+Y(134))
*
*    Y(163) ICaT C0 channel variable
*
      YDOT(163) = k_CaTmv*Y(164) - 4.0*k_CaTpv*Y(163)
     *          + k_CaTmi*Y(169)/h_CaT**3 - k_CaTpi*Y(163)*f_CaT**3
*
*    Y(164) ICaT C1 channel variable
*
      YDOT(164) = 4.0*k_CaTpv*Y(163) - k_CaTmv*Y(164)
     *          + 2.0*k_CaTmv*Y(165) - 3.0*k_CaTpv*Y(164)
     *          + k_CaTmi*Y(170)/h_CaT**2 - k_CaTpi*Y(164)*f_CaT**2
*
*    Y(165) ICaT C2 channel variable
*
      YDOT(165) = 3.0*k_CaTpv*Y(164) - 2.0*k_CaTmv*Y(165)
     *          + 3.0*k_CaTmv*Y(166) - 2.0*k_CaTpv*Y(165)
     *          + k_CaTmi*Y(171)/h_CaT - k_CaTpi*Y(165)*f_CaT
*
*    Y(166) ICaT C3 channel variable
*
      YDOT(166) = 2.0*k_CaTpv*Y(165) - 3.0*k_CaTmv*Y(166)
     *          + 4.0*k_CaTmv*Y(167) - k_CaTpv*Y(166)
     *          + k_CaTmi*Y(172) - k_CaTpi*Y(166)
*
*    Y(167) ICaT C4 channel variable
*
      YDOT(167) = k_CaTpv*Y(166) - 4.0*k_CaTmv*Y(167)
     *          + k_CaTmo*Y(168) - k_CaTpo*Y(167)
     *          + k_CaTmi*Y(173) - k_CaTpi*Y(167)
*
*    Y(168) ICaT O channel variable
*
      YDOT(168) = k_CaTpo*Y(167) - k_CaTmo*Y(168)
     *          + k_CaTmi*Y(174) - k_CaTpi*Y(168)
*
*    Y(169) ICaT I0 channel variable
*
      YDOT(169) = k_CaTmv*Y(170)*h_CaT - 4.0*k_CaTpv*Y(169)/f_CaT
     *          + k_CaTpi*Y(163)*f_CaT**3 - k_CaTmi*Y(169)/h_CaT**3
*
*    Y(170) ICaT I1 channel variable
*
      YDOT(170) = 4.0*k_CaTpv*Y(169)/f_CaT - k_CaTmv*Y(170)*h_CaT
     *          + 2.0*k_CaTmv*Y(171)*h_CaT - 3.0*k_CaTpv*Y(170)/f_CaT
     *          + k_CaTpi*Y(164)*f_CaT**2 - k_CaTmi*Y(170)/h_CaT**2
*
*    Y(171) ICaT I2 channel variable
*
      YDOT(171) = 3.0*k_CaTpv*Y(170)/f_CaT - 2.0*k_CaTmv*Y(171)*h_CaT
     *          + 3.0*k_CaTmv*Y(172)*h_CaT - 2.0*k_CaTpv*Y(171)/f_CaT
     *          + k_CaTpi*Y(165)*f_CaT - k_CaTmi*Y(171)/h_CaT
*
*    Y(172) ICaT I3 channel variable
*
      YDOT(172) = 2.0*k_CaTpv*Y(171)/f_CaT - 3.0*k_CaTmv*Y(172)*h_CaT
     *          + 4.0*k_CaTmv*Y(173) - k_CaTpv*Y(172)
     *          + k_CaTpi*Y(166) - k_CaTmi*Y(172)
*
*     Y(173) ICaT I4 channel variable
*
      YDOT(173) = k_CaTpv*Y(172) - 4.0*k_CaTmv*Y(173)
     *          + k_CaTmo*Y(174) - k_CaTpo*Y(173)
     *          + k_CaTpi*Y(167) - k_CaTmi*Y(173)
*
*    Y(174) ICaT IO channel variable
*
      YDOT(174) = k_CaTpo*Y(173) - k_CaTmo*Y(174)
     *          + k_CaTpi*Y(168) - k_CaTmi*Y(174)
*
*    Y(175) ICaT fraction phophorylated caveolae
*
      YDOT(175) = kIcatPKA*Y(83)*(1-Y(175))/(KIcatmPKA+(1-Y(175)))
     *          - kIcatPP*PPcav*Y(175)/(KIcatmPP+Y(175))
*
*    Y(176) ICaT fraction phophorylated cytosol
*
      YDOT(176) = kIcatPKA*Y(95)*(1-Y(176))/(KIcatmPKA+(1-Y(176)))
     *          - kIcatPP*(PP2Acyt+pp1cytf)*Y(176)/(KIcatmPP+Y(176))
*
*  Small-conductance calcium activated K+ channel
*
      YDOT(177) = k_CaKb1*Y(178) - k_CaKf1*Y(177)
      YDOT(178) = k_CaKf1*Y(177) - k_CaKb1*Y(178)
     *          + k_CaKb2*Y(179) - k_CaKf2*Y(178)
      YDOT(179) = k_CaKf2*Y(178) - k_CaKb2*Y(179)
     *          + k_CaKb3*Y(180) - k_CaKf3*Y(179)
     *          + k_CaKbo1*Y(181) - k_CaKfo1*Y(179)
      YDOT(180) = k_CaKf3*Y(179) - k_CaKb3*Y(180)
     *          + k_CaKbo2*Y(182) - k_CaKfo2*Y(180)
      YDOT(181) = k_CaKfo1*Y(179) - k_CaKbo1*Y(181)
      YDOT(182) = k_CaKfo2*Y(180) - k_CaKbo2*Y(182)
*
*    Y(183) ICaK fraction phophorylated caveolae
*
      YDOT(183) = kIcakPKA*Y(83)*(1-Y(183))/(KIcakmPKA+(1-Y(183)))
     *          - kIcakPP*PPcav*Y(183)/(KIcakmPP+Y(183))
*
  44  format(f12.5,30e15.6)
*
  11  continue
      return
      end
C***********************************************************************
C
      SUBROUTINE RK4(FUN,NEQ,Y,YDOT,TIME,delt,Istim)
C
C***********************************************************************
      implicit none
      integer NEQ,i
      real time,Istim
      real XK1(300),XK2(300),XK3(300),XK4(300),Y1(300),TT,delt
      real Y(300),YDOT(300)
      EXTERNAL FUN
C
          CALL FUN(NEQ,TIME,Y,YDOT,Istim)
C
        DO I=1,NEQ
        XK1(I) = YDOT(I) * delt
        Y1(I) = Y(I) + XK1(I)/2.
        END DO
C
        TT = TIME + delt/2.
          CALL FUN(NEQ,TT,Y1,YDOT,Istim)
C
        DO I=1,NEQ
        XK2(I) = YDOT(I) * delt
        Y1(I) = Y(I) + XK2(I)/2.
        END DO
C
          CALL FUN(NEQ,TT,Y1,YDOT,Istim)
C
        DO I=1,NEQ
        XK3(I) = YDOT(I) * delt
        Y1(I) = Y(I) + XK3(I)
        END DO
C
        TT = TIME + delt
          CALL FUN(NEQ,TT,Y1,YDOT,Istim)
C
        DO I=1,NEQ
        XK4(I) = YDOT(I) * delt
        END DO
C
        DO I=1,NEQ
        Y(I) = Y(I) + (XK1(I)+2.0*XK2(I)+2.0*XK3(I)+XK4(I))/6.0
        END DO
C
        TIME = TIME + delt
        RETURN
        END
C***********************************************************************
