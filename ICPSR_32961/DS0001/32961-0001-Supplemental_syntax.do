/*-------------------------------------------------------------------------*
 |                                                                         
 |            STATA SUPPLEMENTAL SYNTAX FILE FOR ICPSR 32961
 |           STUDY OF WOMEN'S HEALTH ACROSS THE NATION (SWAN),
 |                      2006-2008: VISIT 10 DATASET
 |
 |
 | This Stata missing value recode program is provided for optional use with
 | the Stata system version of this data file as distributed by ICPSR.
 | The program replaces user-defined numeric missing values (e.g., -9)
 | with generic system missing "."  Note that Stata allows you to specify
 | up to 27 unique missing value codes.  Only variables with user-defined
 | missing values are included in this program.
 |
 | To apply the missing value recodes, users need to first open the
 | Stata data file on their system, apply the missing value recodes if
 | desired, then save a new copy of the data file with the missing values
 | applied.  Users are strongly advised to use a different filename when
 | saving the new file.
 |
 *------------------------------------------------------------------------*/

replace LANGINT10 = . if (LANGINT10 <= -1 )
replace PRGNAN10 = . if (PRGNAN10 <= -1 )
replace ANTICO110 = . if (ANTICO110 <= -1 )
replace ACOATW110 = . if (ACOATW110 <= -1 )
replace ANTICO210 = . if (ANTICO210 <= -1 )
replace ACOATW210 = . if (ACOATW210 <= -1 )
replace HEART110 = . if (HEART110 <= -1 )
replace HARTTW110 = . if (HARTTW110 <= -1 )
replace HEART210 = . if (HEART210 <= -1 )
replace HARTTW210 = . if (HARTTW210 <= -1 )
replace CHOLST110 = . if (CHOLST110 <= -1 )
replace CHOLTW110 = . if (CHOLTW110 <= -1 )
replace CHOLST210 = . if (CHOLST210 <= -1 )
replace CHOLTW210 = . if (CHOLTW210 <= -1 )
replace BP110 = . if (BP110 <= -1 )
replace BPTW110 = . if (BPTW110 <= -1 )
replace BP210 = . if (BP210 <= -1 )
replace BPTW210 = . if (BPTW210 <= -1 )
replace DIURET110 = . if (DIURET110 <= -1 )
replace DIURTW110 = . if (DIURTW110 <= -1 )
replace DIURET210 = . if (DIURET210 <= -1 )
replace DIURTW210 = . if (DIURTW210 <= -1 )
replace THYROI110 = . if (THYROI110 <= -1 )
replace THYRTW110 = . if (THYRTW110 <= -1 )
replace THYROI210 = . if (THYROI210 <= -1 )
replace THYRTW210 = . if (THYRTW210 <= -1 )
replace INSULN110 = . if (INSULN110 <= -1 )
replace INSUTW110 = . if (INSUTW110 <= -1 )
replace INSULN210 = . if (INSULN210 <= -1 )
replace INSUTW210 = . if (INSUTW210 <= -1 )
replace NERVS110 = . if (NERVS110 <= -1 )
replace NERVTW110 = . if (NERVTW110 <= -1 )
replace NERVS210 = . if (NERVS210 <= -1 )
replace NERVTW210 = . if (NERVTW210 <= -1 )
replace STEROI110 = . if (STEROI110 <= -1 )
replace STERTW110 = . if (STERTW110 <= -1 )
replace STEROI210 = . if (STEROI210 <= -1 )
replace STERTW210 = . if (STERTW210 <= -1 )
replace ARTHRT110 = . if (ARTHRT110 <= -1 )
replace ARTHTW110 = . if (ARTHTW110 <= -1 )
replace ARTHRT210 = . if (ARTHRT210 <= -1 )
replace ARTHTW210 = . if (ARTHTW210 <= -1 )
replace FERTIL110 = . if (FERTIL110 <= -1 )
replace FRTLTW110 = . if (FRTLTW110 <= -1 )
replace FERTIL210 = . if (FERTIL210 <= -1 )
replace FRTLTW210 = . if (FRTLTW210 <= -1 )
replace BCP110 = . if (BCP110 <= -1 )
replace BCPTWI110 = . if (BCPTWI110 <= -1 )
replace BCP210 = . if (BCP210 <= -1 )
replace BCPTWI210 = . if (BCPTWI210 <= -1 )
replace BCREAS10 = . if (BCREAS10 <= -1 )
replace ESTROG110 = . if (ESTROG110 <= -1 )
replace ESTRTW110 = . if (ESTRTW110 <= -1 )
replace ESTROG210 = . if (ESTROG210 <= -1 )
replace ESTRTW210 = . if (ESTRTW210 <= -1 )
replace ESTRDA110 = . if (ESTRDA110 <= -1 )
replace ESTRDA210 = . if (ESTRDA210 <= -1 )
replace ESTRNJ110 = . if (ESTRNJ110 <= -1 )
replace EINJTW110 = . if (EINJTW110 <= -1 )
replace ESTRNJ210 = . if (ESTRNJ210 <= -1 )
replace EINJTW210 = . if (EINJTW210 <= -1 )
replace COMBIN110 = . if (COMBIN110 <= -1 )
replace COMBTW110 = . if (COMBTW110 <= -1 )
replace COMBIN210 = . if (COMBIN210 <= -1 )
replace COMBTW210 = . if (COMBTW210 <= -1 )
replace PROGES110 = . if (PROGES110 <= -1 )
replace PROGTW110 = . if (PROGTW110 <= -1 )
replace PROGES210 = . if (PROGES210 <= -1 )
replace PROGTW210 = . if (PROGTW210 <= -1 )
replace PROGDA110 = . if (PROGDA110 <= -1 )
replace PROGDA210 = . if (PROGDA210 <= -1 )
replace OHRM_110 = . if (OHRM_110 <= -1 )
replace OHRMMO110 = . if (OHRMMO110 <= -1 )
replace OHRM_210 = . if (OHRM_210 <= -1 )
replace OHRMMO210 = . if (OHRMMO210 <= -1 )
replace OHRM_310 = . if (OHRM_310 <= -1 )
replace OHRMMO310 = . if (OHRMMO310 <= -1 )
replace OHRM_410 = . if (OHRM_410 <= -1 )
replace OHRMMO410 = . if (OHRMMO410 <= -1 )
replace OSTEIV110 = . if (OSTEIV110 <= -1 )
replace OSTIVL110 = . if (OSTIVL110 <= -1 )
replace OSTEON110 = . if (OSTEON110 <= -1 )
replace OSTNTW110 = . if (OSTNTW110 <= -1 )
replace OSTN1W110 = . if (OSTN1W110 <= -1 )
replace OSTN1MO110 = . if (OSTN1MO110 <= -1 )
replace OSTEON210 = . if (OSTEON210 <= -1 )
replace OSTNTW210 = . if (OSTNTW210 <= -1 )
replace OSTN1W210 = . if (OSTN1W210 <= -1 )
replace OSTN1MO210 = . if (OSTN1MO210 <= -1 )
replace OTHMED110 = . if (OTHMED110 <= -1 )
replace OTHRTW110 = . if (OTHRTW110 <= -1 )
replace OTHMED210 = . if (OTHMED210 <= -1 )
replace OTHRTW210 = . if (OTHRTW210 <= -1 )
replace OTHMED310 = . if (OTHMED310 <= -1 )
replace OTHRTW310 = . if (OTHRTW310 <= -1 )
replace OTHMED410 = . if (OTHMED410 <= -1 )
replace OTHRTW410 = . if (OTHRTW410 <= -1 )
replace OTHMED510 = . if (OTHMED510 <= -1 )
replace OTHRTW510 = . if (OTHRTW510 <= -1 )
replace OTHMED610 = . if (OTHMED610 <= -1 )
replace OTHRTW610 = . if (OTHRTW610 <= -1 )
replace OTHMED710 = . if (OTHMED710 <= -1 )
replace OTHRTW710 = . if (OTHRTW710 <= -1 )
replace OTHMED810 = . if (OTHMED810 <= -1 )
replace OTHRTW810 = . if (OTHRTW810 <= -1 )
replace OTHMED910 = . if (OTHMED910 <= -1 )
replace OTHRTW910 = . if (OTHRTW910 <= -1 )
replace OTHME1010 = . if (OTHME1010 <= -1 )
replace OTHTW1010 = . if (OTHTW1010 <= -1 )
replace OTHME1110 = . if (OTHME1110 <= -1 )
replace OTHTW1110 = . if (OTHTW1110 <= -1 )
replace OTHME1210 = . if (OTHME1210 <= -1 )
replace OTHTW1210 = . if (OTHTW1210 <= -1 )
replace OTHME1310 = . if (OTHME1310 <= -1 )
replace OTHTW1310 = . if (OTHTW1310 <= -1 )
replace OTHME1410 = . if (OTHME1410 <= -1 )
replace OTHTW1410 = . if (OTHTW1410 <= -1 )
replace OTHME1510 = . if (OTHME1510 <= -1 )
replace OTHTW1510 = . if (OTHTW1510 <= -1 )
replace ESTLSTV10 = . if (ESTLSTV10 <= -1 )
replace REDUHAR10 = . if (REDUHAR10 <= -1 )
replace OSTEOPO10 = . if (OSTEOPO10 <= -1 )
replace MENOSYM10 = . if (MENOSYM10 <= -1 )
replace YOUNGLK10 = . if (YOUNGLK10 <= -1 )
replace HCPADVI10 = . if (HCPADVI10 <= -1 )
replace FRNADVI10 = . if (FRNADVI10 <= -1 )
replace IMPRMEM10 = . if (IMPRMEM10 <= -1 )
replace REGPERI10 = . if (REGPERI10 <= -1 )
replace HORMOTH10 = . if (HORMOTH10 <= -1 )
replace DONTKNO10 = . if (DONTKNO10 <= -1 )
replace PRBBLEE10 = . if (PRBBLEE10 <= -1 )
replace HAVEPER10 = . if (HAVEPER10 <= -1 )
replace LIKEFEL10 = . if (LIKEFEL10 <= -1 )
replace SIDEEFF10 = . if (SIDEEFF10 <= -1 )
replace CANCER10 = . if (CANCER10 <= -1 )
replace ADVISTO10 = . if (ADVISTO10 <= -1 )
replace EXPENSI10 = . if (EXPENSI10 <= -1 )
replace NOLIKE10 = . if (NOLIKE10 <= -1 )
replace NOREMEB10 = . if (NOREMEB10 <= -1 )
replace DNTKNOW10 = . if (DNTKNOW10 <= -1 )
replace STOPOTH10 = . if (STOPOTH10 <= -1 )
replace NOREASO10 = . if (NOREASO10 <= -1 )
replace NEWSRPT10 = . if (NEWSRPT10 <= -1 )
replace REGVITA10 = . if (REGVITA10 <= -1 )
replace ONCEADA10 = . if (ONCEADA10 <= -1 )
replace ANTIOXI10 = . if (ANTIOXI10 <= -1 )
replace VITCOMB10 = . if (VITCOMB10 <= -1 )
replace VTMOTH110 = . if (VTMOTH110 <= -1 )
replace VTMOTH210 = . if (VTMOTH210 <= -1 )
replace VTMOTH310 = . if (VTMOTH310 <= -1 )
replace VTMOTH410 = . if (VTMOTH410 <= -1 )
replace VITAMNA10 = . if (VITAMNA10 <= -1 )
replace BETACAR10 = . if (BETACAR10 <= -1 )
replace VITAMNC10 = . if (VITAMNC10 <= -1 )
replace VITAMND10 = . if (VITAMND10 <= -1 )
replace VITAMNE10 = . if (VITAMNE10 <= -1 )
replace CALCTUM10 = . if (CALCTUM10 <= -1 )
replace IRON10 = . if (IRON10 <= -1 )
replace ZINC10 = . if (ZINC10 <= -1 )
replace SELENIU10 = . if (SELENIU10 <= -1 )
replace FOLATE10 = . if (FOLATE10 <= -1 )
replace VTMSING10 = . if (VTMSING10 <= -1 )
replace SVTMOT110 = . if (SVTMOT110 <= -1 )
replace SVTMOT210 = . if (SVTMOT210 <= -1 )
replace SVTMOT310 = . if (SVTMOT310 <= -1 )
replace SVTMOT410 = . if (SVTMOT410 <= -1 )
replace SVTMOT510 = . if (SVTMOT510 <= -1 )
replace SVTMOT610 = . if (SVTMOT610 <= -1 )
replace SVTMOT710 = . if (SVTMOT710 <= -1 )
replace SVTMOT810 = . if (SVTMOT810 <= -1 )
replace SVTMOT910 = . if (SVTMOT910 <= -1 )
replace SVTMO1010 = . if (SVTMO1010 <= -1 )
replace PAIN110 = . if (PAIN110 <= -1 )
replace PAINTW110 = . if (PAINTW110 <= -1 )
replace PAIN210 = . if (PAIN210 <= -1 )
replace PAINTW210 = . if (PAINTW210 <= -1 )
replace SLEEP110 = . if (SLEEP110 <= -1 )
replace SLEPTW110 = . if (SLEPTW110 <= -1 )
replace SLEEP210 = . if (SLEEP210 <= -1 )
replace SLEPTW210 = . if (SLEPTW210 <= -1 )
replace OTC110 = . if (OTC110 <= -1 )
replace OTCTW110 = . if (OTCTW110 <= -1 )
replace OTC210 = . if (OTC210 <= -1 )
replace OTCTW210 = . if (OTCTW210 <= -1 )
replace OTC310 = . if (OTC310 <= -1 )
replace OTCTW310 = . if (OTCTW310 <= -1 )
replace OTC410 = . if (OTC410 <= -1 )
replace OTCTW410 = . if (OTCTW410 <= -1 )
replace OTC510 = . if (OTC510 <= -1 )
replace OTCTW510 = . if (OTCTW510 <= -1 )
replace OTC610 = . if (OTC610 <= -1 )
replace OTCTW610 = . if (OTCTW610 <= -1 )
replace OTC710 = . if (OTC710 <= -1 )
replace OTCTW710 = . if (OTCTW710 <= -1 )
replace OTC810 = . if (OTC810 <= -1 )
replace OTCTW810 = . if (OTCTW810 <= -1 )
replace OTC910 = . if (OTC910 <= -1 )
replace OTCTW910 = . if (OTCTW910 <= -1 )
replace OTC1010 = . if (OTC1010 <= -1 )
replace OTCTW1010 = . if (OTCTW1010 <= -1 )
replace OTC1110 = . if (OTC1110 <= -1 )
replace OTCTW1110 = . if (OTCTW1110 <= -1 )
replace OTC1210 = . if (OTC1210 <= -1 )
replace OTCTW1210 = . if (OTCTW1210 <= -1 )
replace OTC1310 = . if (OTC1310 <= -1 )
replace OTCTW1310 = . if (OTCTW1310 <= -1 )
replace OTC1410 = . if (OTC1410 <= -1 )
replace OTCTW1410 = . if (OTCTW1410 <= -1 )
replace OTC1510 = . if (OTC1510 <= -1 )
replace OTCTW1510 = . if (OTCTW1510 <= -1 )
replace SOYYSNO10 = . if (SOYYSNO10 <= -1 )
replace SOYPROT10 = . if (SOYPROT10 <= -1 )
replace CEREACA10 = . if (CEREACA10 <= -1 )
replace BREADCA10 = . if (BREADCA10 <= -1 )
replace ORANGCA10 = . if (ORANGCA10 <= -1 )
replace ACUPUNC10 = . if (ACUPUNC10 <= -1 )
replace ACUPHAR10 = . if (ACUPHAR10 <= -1 )
replace ACUPOST10 = . if (ACUPOST10 <= -1 )
replace ACUPMEN10 = . if (ACUPMEN10 <= -1 )
replace ACUPLOO10 = . if (ACUPLOO10 <= -1 )
replace ACUPMEM10 = . if (ACUPMEM10 <= -1 )
replace ACUPPER10 = . if (ACUPPER10 <= -1 )
replace ACUPGEN10 = . if (ACUPGEN10 <= -1 )
replace ACUPWGH10 = . if (ACUPWGH10 <= -1 )
replace ACUPADV10 = . if (ACUPADV10 <= -1 )
replace ACUPOTH10 = . if (ACUPOTH10 <= -1 )
replace BCOHOSH10 = . if (BCOHOSH10 <= -1 )
replace BCOHHAR10 = . if (BCOHHAR10 <= -1 )
replace BCOHOST10 = . if (BCOHOST10 <= -1 )
replace BCOHMEN10 = . if (BCOHMEN10 <= -1 )
replace BCOHLOO10 = . if (BCOHLOO10 <= -1 )
replace BCOHMEM10 = . if (BCOHMEM10 <= -1 )
replace BCOHPER10 = . if (BCOHPER10 <= -1 )
replace BCOHGEN10 = . if (BCOHGEN10 <= -1 )
replace BCOHWGH10 = . if (BCOHWGH10 <= -1 )
replace BCOHADV10 = . if (BCOHADV10 <= -1 )
replace BCOHOTH10 = . if (BCOHOTH10 <= -1 )
replace DQUAI10 = . if (DQUAI10 <= -1 )
replace DQUAHAR10 = . if (DQUAHAR10 <= -1 )
replace DQUAOST10 = . if (DQUAOST10 <= -1 )
replace DQUAMEN10 = . if (DQUAMEN10 <= -1 )
replace DQUALOO10 = . if (DQUALOO10 <= -1 )
replace DQUAMEM10 = . if (DQUAMEM10 <= -1 )
replace DQUAPER10 = . if (DQUAPER10 <= -1 )
replace DQUAGEN10 = . if (DQUAGEN10 <= -1 )
replace DQUAWGH10 = . if (DQUAWGH10 <= -1 )
replace DQUAADV10 = . if (DQUAADV10 <= -1 )
replace DQUAOTH10 = . if (DQUAOTH10 <= -1 )
replace DIETNUT10 = . if (DIETNUT10 <= -1 )
replace DIETHAR10 = . if (DIETHAR10 <= -1 )
replace DIETOST10 = . if (DIETOST10 <= -1 )
replace DIETMEN10 = . if (DIETMEN10 <= -1 )
replace DIETLOO10 = . if (DIETLOO10 <= -1 )
replace DIETMEM10 = . if (DIETMEM10 <= -1 )
replace DIETPER10 = . if (DIETPER10 <= -1 )
replace DIETGEN10 = . if (DIETGEN10 <= -1 )
replace DIETWGH10 = . if (DIETWGH10 <= -1 )
replace DIETADV10 = . if (DIETADV10 <= -1 )
replace DIETOTH10 = . if (DIETOTH10 <= -1 )
replace EXERCIS10 = . if (EXERCIS10 <= -1 )
replace EXERHAR10 = . if (EXERHAR10 <= -1 )
replace EXEROST10 = . if (EXEROST10 <= -1 )
replace EXERMEN10 = . if (EXERMEN10 <= -1 )
replace EXERLOO10 = . if (EXERLOO10 <= -1 )
replace EXERMEM10 = . if (EXERMEM10 <= -1 )
replace EXERPER10 = . if (EXERPER10 <= -1 )
replace EXERGEN10 = . if (EXERGEN10 <= -1 )
replace EXERWGH10 = . if (EXERWGH10 <= -1 )
replace EXERADV10 = . if (EXERADV10 <= -1 )
replace EXEROTH10 = . if (EXEROTH10 <= -1 )
replace FLAXSEE10 = . if (FLAXSEE10 <= -1 )
replace FLAXHAR10 = . if (FLAXHAR10 <= -1 )
replace FLAXOST10 = . if (FLAXOST10 <= -1 )
replace FLAXMEN10 = . if (FLAXMEN10 <= -1 )
replace FLAXLOO10 = . if (FLAXLOO10 <= -1 )
replace FLAXMEM10 = . if (FLAXMEM10 <= -1 )
replace FLAXPER10 = . if (FLAXPER10 <= -1 )
replace FLAXGEN10 = . if (FLAXGEN10 <= -1 )
replace FLAXWGH10 = . if (FLAXWGH10 <= -1 )
replace FLAXADV10 = . if (FLAXADV10 <= -1 )
replace FLAXOTH10 = . if (FLAXOTH10 <= -1 )
replace GINKGO10 = . if (GINKGO10 <= -1 )
replace GINKHAR10 = . if (GINKHAR10 <= -1 )
replace GINKOST10 = . if (GINKOST10 <= -1 )
replace GINKMEN10 = . if (GINKMEN10 <= -1 )
replace GINKLOO10 = . if (GINKLOO10 <= -1 )
replace GINKMEM10 = . if (GINKMEM10 <= -1 )
replace GINKPER10 = . if (GINKPER10 <= -1 )
replace GINKGEN10 = . if (GINKGEN10 <= -1 )
replace GINKWGH10 = . if (GINKWGH10 <= -1 )
replace GINKADV10 = . if (GINKADV10 <= -1 )
replace GINKOTH10 = . if (GINKOTH10 <= -1 )
replace GINSENG10 = . if (GINSENG10 <= -1 )
replace GINSHAR10 = . if (GINSHAR10 <= -1 )
replace GINSOST10 = . if (GINSOST10 <= -1 )
replace GINSMEN10 = . if (GINSMEN10 <= -1 )
replace GINSLOO10 = . if (GINSLOO10 <= -1 )
replace GINSMEM10 = . if (GINSMEM10 <= -1 )
replace GINSPER10 = . if (GINSPER10 <= -1 )
replace GINSGEN10 = . if (GINSGEN10 <= -1 )
replace GINSWGH10 = . if (GINSWGH10 <= -1 )
replace GINSADV10 = . if (GINSADV10 <= -1 )
replace GINSOTH10 = . if (GINSOTH10 <= -1 )
replace GLUSAMI10 = . if (GLUSAMI10 <= -1 )
replace GLUSHAR10 = . if (GLUSHAR10 <= -1 )
replace GLUSOST10 = . if (GLUSOST10 <= -1 )
replace GLUSMEN10 = . if (GLUSMEN10 <= -1 )
replace GLUSLOO10 = . if (GLUSLOO10 <= -1 )
replace GLUSMEM10 = . if (GLUSMEM10 <= -1 )
replace GLUSPER10 = . if (GLUSPER10 <= -1 )
replace GLUSGEN10 = . if (GLUSGEN10 <= -1 )
replace GLUSWGH10 = . if (GLUSWGH10 <= -1 )
replace GLUSADV10 = . if (GLUSADV10 <= -1 )
replace GLUSOTH10 = . if (GLUSOTH10 <= -1 )
replace MYAMPRO10 = . if (MYAMPRO10 <= -1 )
replace MYAMHAR10 = . if (MYAMHAR10 <= -1 )
replace MYAMOST10 = . if (MYAMOST10 <= -1 )
replace MYAMMEN10 = . if (MYAMMEN10 <= -1 )
replace MYAMLOO10 = . if (MYAMLOO10 <= -1 )
replace MYAMMEM10 = . if (MYAMMEM10 <= -1 )
replace MYAMPER10 = . if (MYAMPER10 <= -1 )
replace MYAMGEN10 = . if (MYAMGEN10 <= -1 )
replace MYAMWGH10 = . if (MYAMWGH10 <= -1 )
replace MYAMADV10 = . if (MYAMADV10 <= -1 )
replace MYAMOTH10 = . if (MYAMOTH10 <= -1 )
replace PRAYER10 = . if (PRAYER10 <= -1 )
replace PRAYHAR10 = . if (PRAYHAR10 <= -1 )
replace PRAYOST10 = . if (PRAYOST10 <= -1 )
replace PRAYMEN10 = . if (PRAYMEN10 <= -1 )
replace PRAYLOO10 = . if (PRAYLOO10 <= -1 )
replace PRAYMEM10 = . if (PRAYMEM10 <= -1 )
replace PRAYPER10 = . if (PRAYPER10 <= -1 )
replace PRAYGEN10 = . if (PRAYGEN10 <= -1 )
replace PRAYWGH10 = . if (PRAYWGH10 <= -1 )
replace PRAYADV10 = . if (PRAYADV10 <= -1 )
replace PRAYOTH10 = . if (PRAYOTH10 <= -1 )
replace SELFHEL10 = . if (SELFHEL10 <= -1 )
replace SELFHAR10 = . if (SELFHAR10 <= -1 )
replace SELFOST10 = . if (SELFOST10 <= -1 )
replace SELFMEN10 = . if (SELFMEN10 <= -1 )
replace SELFLOO10 = . if (SELFLOO10 <= -1 )
replace SELFMEM10 = . if (SELFMEM10 <= -1 )
replace SELFPER10 = . if (SELFPER10 <= -1 )
replace SELFGEN10 = . if (SELFGEN10 <= -1 )
replace SELFWGH10 = . if (SELFWGH10 <= -1 )
replace SELFADV10 = . if (SELFADV10 <= -1 )
replace SELFOTH10 = . if (SELFOTH10 <= -1 )
replace SOYSUPP10 = . if (SOYSUPP10 <= -1 )
replace SOYHAR10 = . if (SOYHAR10 <= -1 )
replace SOYOST10 = . if (SOYOST10 <= -1 )
replace SOYMEN10 = . if (SOYMEN10 <= -1 )
replace SOYLOO10 = . if (SOYLOO10 <= -1 )
replace SOYMEM10 = . if (SOYMEM10 <= -1 )
replace SOYPER10 = . if (SOYPER10 <= -1 )
replace SOYGEN10 = . if (SOYGEN10 <= -1 )
replace SOYWGH10 = . if (SOYWGH10 <= -1 )
replace SOYADV10 = . if (SOYADV10 <= -1 )
replace SOYOTH10 = . if (SOYOTH10 <= -1 )
replace WORTSTJ10 = . if (WORTSTJ10 <= -1 )
replace WORTHAR10 = . if (WORTHAR10 <= -1 )
replace WORTOST10 = . if (WORTOST10 <= -1 )
replace WORTMEN10 = . if (WORTMEN10 <= -1 )
replace WORTLOO10 = . if (WORTLOO10 <= -1 )
replace WORTMEM10 = . if (WORTMEM10 <= -1 )
replace WORTPER10 = . if (WORTPER10 <= -1 )
replace WORTGEN10 = . if (WORTGEN10 <= -1 )
replace WORTWGH10 = . if (WORTWGH10 <= -1 )
replace WORTADV10 = . if (WORTADV10 <= -1 )
replace WORTOTH10 = . if (WORTOTH10 <= -1 )
replace WVITAMI10 = . if (WVITAMI10 <= -1 )
replace WVITHAR10 = . if (WVITHAR10 <= -1 )
replace WVITOST10 = . if (WVITOST10 <= -1 )
replace WVITMEN10 = . if (WVITMEN10 <= -1 )
replace WVITLOO10 = . if (WVITLOO10 <= -1 )
replace WVITMEM10 = . if (WVITMEM10 <= -1 )
replace WVITPER10 = . if (WVITPER10 <= -1 )
replace WVITGEN10 = . if (WVITGEN10 <= -1 )
replace WVITWGH10 = . if (WVITWGH10 <= -1 )
replace WVITADV10 = . if (WVITADV10 <= -1 )
replace WVITOTH10 = . if (WVITOTH10 <= -1 )
replace YOGA10 = . if (YOGA10 <= -1 )
replace YOGAHAR10 = . if (YOGAHAR10 <= -1 )
replace YOGAOST10 = . if (YOGAOST10 <= -1 )
replace YOGAMEN10 = . if (YOGAMEN10 <= -1 )
replace YOGALOO10 = . if (YOGALOO10 <= -1 )
replace YOGAMEM10 = . if (YOGAMEM10 <= -1 )
replace YOGAPER10 = . if (YOGAPER10 <= -1 )
replace YOGAGEN10 = . if (YOGAGEN10 <= -1 )
replace YOGAWGH10 = . if (YOGAWGH10 <= -1 )
replace YOGAADV10 = . if (YOGAADV10 <= -1 )
replace YOGAOTH10 = . if (YOGAOTH10 <= -1 )
replace BOTANIC10 = . if (BOTANIC10 <= -1 )
replace BOTAHAR10 = . if (BOTAHAR10 <= -1 )
replace BOTAOST10 = . if (BOTAOST10 <= -1 )
replace BOTAMEN10 = . if (BOTAMEN10 <= -1 )
replace BOTALOO10 = . if (BOTALOO10 <= -1 )
replace BOTAMEM10 = . if (BOTAMEM10 <= -1 )
replace BOTAPER10 = . if (BOTAPER10 <= -1 )
replace BOTAGEN10 = . if (BOTAGEN10 <= -1 )
replace BOTAWGH10 = . if (BOTAWGH10 <= -1 )
replace BOTAADV10 = . if (BOTAADV10 <= -1 )
replace BOTAOTH10 = . if (BOTAOTH10 <= -1 )
replace HERBALT10 = . if (HERBALT10 <= -1 )
replace HTEAHAR10 = . if (HTEAHAR10 <= -1 )
replace HTEAOST10 = . if (HTEAOST10 <= -1 )
replace HTEAMEN10 = . if (HTEAMEN10 <= -1 )
replace HTEALOO10 = . if (HTEALOO10 <= -1 )
replace HTEAMEM10 = . if (HTEAMEM10 <= -1 )
replace HTEAPER10 = . if (HTEAPER10 <= -1 )
replace HTEAGEN10 = . if (HTEAGEN10 <= -1 )
replace HTEAWGH10 = . if (HTEAWGH10 <= -1 )
replace HTEAADV10 = . if (HTEAADV10 <= -1 )
replace HTEAOTH10 = . if (HTEAOTH10 <= -1 )
replace OTHALT10 = . if (OTHALT10 <= -1 )
replace OTHHAR10 = . if (OTHHAR10 <= -1 )
replace OTHOST10 = . if (OTHOST10 <= -1 )
replace OTHMEN10 = . if (OTHMEN10 <= -1 )
replace OTHLOO10 = . if (OTHLOO10 <= -1 )
replace OTHMEM10 = . if (OTHMEM10 <= -1 )
replace OTHPER10 = . if (OTHPER10 <= -1 )
replace OTHGEN10 = . if (OTHGEN10 <= -1 )
replace OTHWGH10 = . if (OTHWGH10 <= -1 )
replace OTHADV10 = . if (OTHADV10 <= -1 )
replace OTHALTR10 = . if (OTHALTR10 <= -1 )
replace OTHALT210 = . if (OTHALT210 <= -1 )
replace OT2HAR10 = . if (OT2HAR10 <= -1 )
replace OT2OST10 = . if (OT2OST10 <= -1 )
replace OT2MEN10 = . if (OT2MEN10 <= -1 )
replace OT2LOO10 = . if (OT2LOO10 <= -1 )
replace OT2MEM10 = . if (OT2MEM10 <= -1 )
replace OT2PER10 = . if (OT2PER10 <= -1 )
replace OT2GEN10 = . if (OT2GEN10 <= -1 )
replace OT2WGH10 = . if (OT2WGH10 <= -1 )
replace OT2ADV10 = . if (OT2ADV10 <= -1 )
replace OT2ALT10 = . if (OT2ALT10 <= -1 )
replace OTHALT310 = . if (OTHALT310 <= -1 )
replace OT3HAR10 = . if (OT3HAR10 <= -1 )
replace OT3OST10 = . if (OT3OST10 <= -1 )
replace OT3MEN10 = . if (OT3MEN10 <= -1 )
replace OT3LOO10 = . if (OT3LOO10 <= -1 )
replace OT3MEM10 = . if (OT3MEM10 <= -1 )
replace OT3PER10 = . if (OT3PER10 <= -1 )
replace OT3GEN10 = . if (OT3GEN10 <= -1 )
replace OT3WGH10 = . if (OT3WGH10 <= -1 )
replace OT3ADV10 = . if (OT3ADV10 <= -1 )
replace OT3ALT10 = . if (OT3ALT10 <= -1 )
replace ANEMIA10 = . if (ANEMIA10 <= -1 )
replace DIABETE10 = . if (DIABETE10 <= -1 )
replace HIGHBP10 = . if (HIGHBP10 <= -1 )
replace HBCHOLE10 = . if (HBCHOLE10 <= -1 )
replace MIGRAIN10 = . if (MIGRAIN10 <= -1 )
replace STROKE10 = . if (STROKE10 <= -1 )
replace OSTEOAR10 = . if (OSTEOAR10 <= -1 )
replace THYROID10 = . if (THYROID10 <= -1 )
replace HEARTAT10 = . if (HEARTAT10 <= -1 )
replace ANGINA10 = . if (ANGINA10 <= -1 )
replace OSTEOPR10 = . if (OSTEOPR10 <= -1 )
replace SKCNCER10 = . if (SKCNCER10 <= -1 )
replace MECNCER10 = . if (MECNCER10 <= -1 )
replace NMECNCR10 = . if (NMECNCR10 <= -1 )
replace CANCERS10 = . if (CANCERS10 <= -1 )
replace PSITECA10 = . if (PSITECA10 <= -1 )
replace TAMOXIF10 = . if (TAMOXIF10 <= -1 )
replace CHEMOTH10 = . if (CHEMOTH10 <= -1 )
replace BROKEBO10 = . if (BROKEBO10 <= -1 )
replace HAPPEN110 = . if (HAPPEN110 <= -1 )
replace HAPPEN210 = . if (HAPPEN210 <= -1 )
replace HAPPEN310 = . if (HAPPEN310 <= -1 )
replace DANDC10 = . if (DANDC10 <= -1 )
replace NUMDAND10 = . if (NUMDAND10 <= -1 )
replace HYSTERE10 = . if (HYSTERE10 <= -1 )
replace OOPHORE10 = . if (OOPHORE10 <= -1 )
replace ONEOVAR10 = . if (ONEOVAR10 <= -1 )
replace ABLATIN10 = . if (ABLATIN10 <= -1 )
replace UTERPRO10 = . if (UTERPRO10 <= -1 )
replace THYRREM10 = . if (THYRREM10 <= -1 )
replace ENDO10 = . if (ENDO10 <= -1 )
replace PELVCPN10 = . if (PELVCPN10 <= -1 )
replace PROLAPS10 = . if (PROLAPS10 <= -1 )
replace ABBLEED10 = . if (ABBLEED10 <= -1 )
replace FIBRUTR10 = . if (FIBRUTR10 <= -1 )
replace PRVIDER10 = . if (PRVIDER10 <= -1 )
replace PROFDEG10 = . if (PROFDEG10 <= -1 )
replace PROVSPC10 = . if (PROVSPC10 <= -1 )
replace PROVTIM10 = . if (PROVTIM10 <= -1 )
replace BLEEDNG10 = . if (BLEEDNG10 <= -1 )
replace BLD3MON10 = . if (BLD3MON10 <= -1 )
replace DESCPER10 = . if (DESCPER10 <= -1 )
replace LENGCYL10 = . if (LENGCYL10 <= -1 )
replace PRGNANT10 = . if (PRGNANT10 <= -1 )
replace OUTCOME10 = . if (OUTCOME10 <= -1 )
replace BRSTFEE10 = . if (BRSTFEE10 <= -1 )
replace QLTYLIF10 = . if (QLTYLIF10 <= -1 )
replace LISTEN10 = . if (LISTEN10 <= -1 )
replace TAKETOM10 = . if (TAKETOM10 <= -1 )
replace CONFIDE10 = . if (CONFIDE10 <= -1 )
replace HELPSIC10 = . if (HELPSIC10 <= -1 )
replace CONTROL10 = . if (CONTROL10 <= -1 )
replace ABILITY10 = . if (ABILITY10 <= -1 )
replace YOURWAY10 = . if (YOURWAY10 <= -1 )
replace PILING10 = . if (PILING10 <= -1 )
replace BOTHER10 = . if (BOTHER10 <= -1 )
replace APPETIT10 = . if (APPETIT10 <= -1 )
replace BLUES10 = . if (BLUES10 <= -1 )
replace GOOD10 = . if (GOOD10 <= -1 )
replace KEEPMIN10 = . if (KEEPMIN10 <= -1 )
replace DEPRESS10 = . if (DEPRESS10 <= -1 )
replace EFFORT10 = . if (EFFORT10 <= -1 )
replace HOPEFUL10 = . if (HOPEFUL10 <= -1 )
replace FAILURE10 = . if (FAILURE10 <= -1 )
replace FEARFUL10 = . if (FEARFUL10 <= -1 )
replace RESTLES10 = . if (RESTLES10 <= -1 )
replace HAPPY10 = . if (HAPPY10 <= -1 )
replace TALKLES10 = . if (TALKLES10 <= -1 )
replace LONELY10 = . if (LONELY10 <= -1 )
replace UNFRNDL10 = . if (UNFRNDL10 <= -1 )
replace ENJOY10 = . if (ENJOY10 <= -1 )
replace CRYING10 = . if (CRYING10 <= -1 )
replace SAD10 = . if (SAD10 <= -1 )
replace DISLIKE10 = . if (DISLIKE10 <= -1 )
replace GETGOIN10 = . if (GETGOIN10 <= -1 )
replace CHNGJOB10 = . if (CHNGJOB10 <= -1 )
replace JOB10 = . if (JOB10 <= -1 )
replace CHANGHR10 = . if (CHANGHR10 <= -1 )
replace HOURSPA10 = . if (HOURSPA10 <= -1 )
replace MARITAL10 = . if (MARITAL10 <= -1 )
replace MAILLOC10 = . if (MAILLOC10 <= -1 )
replace BOTHOVR10 = . if (BOTHOVR10 <= -1 )
replace LANGSAA10 = . if (LANGSAA10 <= -1 )
replace OVERHLT10 = . if (OVERHLT10 <= -1 )
replace HOSPSTA10 = . if (HOSPSTA10 <= -1 )
replace MDTALK10 = . if (MDTALK10 <= -1 )
replace NERVES10 = . if (NERVES10 <= -1 )
replace PAPSMEA10 = . if (PAPSMEA10 <= -1 )
replace BRSTEXA10 = . if (BRSTEXA10 <= -1 )
replace MAMOGRA10 = . if (MAMOGRA10 <= -1 )
replace HLTHSER10 = . if (HLTHSER10 <= -1 )
replace INSURAN10 = . if (INSURAN10 <= -1 )
replace NOTAFFR10 = . if (NOTAFFR10 <= -1 )
replace NOTRANS10 = . if (NOTRANS10 <= -1 )
replace NOPROVI10 = . if (NOPROVI10 <= -1 )
replace TOOBUSY10 = . if (TOOBUSY10 <= -1 )
replace NOTRUST10 = . if (NOTRUST10 <= -1 )
replace BETTROF10 = . if (BETTROF10 <= -1 )
replace FAILOTH10 = . if (FAILOTH10 <= -1 )
replace SMOKERE10 = . if (SMOKERE10 <= -1 )
replace AVCIGDA10 = . if (AVCIGDA10 <= -1 )
replace HHMEMSM10 = . if (HHMEMSM10 <= -1 )
replace HOMEXPD10 = . if (HOMEXPD10 <= -1 )
replace HOMEXPH10 = . if (HOMEXPH10 <= -1 )
replace DRNKBEE10 = . if (DRNKBEE10 <= -1 )
replace GLASBEE10 = . if (GLASBEE10 <= -1 )
replace GLASWIN10 = . if (GLASWIN10 <= -1 )
replace GLASLIQ10 = . if (GLASLIQ10 <= -1 )
replace HLTHAYR10 = . if (HLTHAYR10 <= -1 )
replace V_ACTI10 = . if (V_ACTI10 <= -1 )
replace M_ACTI10 = . if (M_ACTI10 <= -1 )
replace LIFTING10 = . if (LIFTING10 <= -1 )
replace CLIMBS10 = . if (CLIMBS10 <= -1 )
replace CLIMB1_10 = . if (CLIMB1_10 <= -1 )
replace BENDING10 = . if (BENDING10 <= -1 )
replace WALKM10 = . if (WALKM10 <= -1 )
replace WALKS10 = . if (WALKS10 <= -1 )
replace WALK1_10 = . if (WALK1_10 <= -1 )
replace BATHING10 = . if (BATHING10 <= -1 )
replace PHYCTDW10 = . if (PHYCTDW10 <= -1 )
replace PHYACCO10 = . if (PHYACCO10 <= -1 )
replace PHYLIMI10 = . if (PHYLIMI10 <= -1 )
replace PHYDFCL10 = . if (PHYDFCL10 <= -1 )
replace EMOCTDW10 = . if (EMOCTDW10 <= -1 )
replace EMOACCO10 = . if (EMOACCO10 <= -1 )
replace EMOCARE10 = . if (EMOCARE10 <= -1 )
replace INTERFR10 = . if (INTERFR10 <= -1 )
replace BODYPAI10 = . if (BODYPAI10 <= -1 )
replace PAINTRF10 = . if (PAINTRF10 <= -1 )
replace PEP10 = . if (PEP10 <= -1 )
replace NERV4WK10 = . if (NERV4WK10 <= -1 )
replace CHER4WK10 = . if (CHER4WK10 <= -1 )
replace CALM4WK10 = . if (CALM4WK10 <= -1 )
replace ENERGY10 = . if (ENERGY10 <= -1 )
replace BLUE4WK10 = . if (BLUE4WK10 <= -1 )
replace WORNOUT10 = . if (WORNOUT10 <= -1 )
replace HAPY4WK10 = . if (HAPY4WK10 <= -1 )
replace TIRED10 = . if (TIRED10 <= -1 )
replace SOCIAL10 = . if (SOCIAL10 <= -1 )
replace HEALSIC10 = . if (HEALSIC10 <= -1 )
replace HEALTHY10 = . if (HEALTHY10 <= -1 )
replace HEALWOR10 = . if (HEALWOR10 <= -1 )
replace HEALEXC10 = . if (HEALEXC10 <= -1 )
replace HOTFLAS10 = . if (HOTFLAS10 <= -1 )
replace NUMHOTF10 = . if (NUMHOTF10 <= -1.0 )
replace BOTHOTF10 = . if (BOTHOTF10 <= -1 )
replace COLDSWE10 = . if (COLDSWE10 <= -1 )
replace NUMCLDS10 = . if (NUMCLDS10 <= -1.0 )
replace BOTCLDS10 = . if (BOTCLDS10 <= -1 )
replace NITESWE10 = . if (NITESWE10 <= -1 )
replace NUMNITS10 = . if (NUMNITS10 <= -1.0 )
replace BOTNITS10 = . if (BOTNITS10 <= -1 )
replace STIFF10 = . if (STIFF10 <= -1 )
replace BOTSTIF10 = . if (BOTSTIF10 <= -1 )
replace IRRITAB10 = . if (IRRITAB10 <= -1 )
replace BOTIRRT10 = . if (BOTIRRT10 <= -1 )
replace NRVOUS10 = . if (NRVOUS10 <= -1 )
replace BOTNERV10 = . if (BOTNERV10 <= -1 )
replace ACHES10 = . if (ACHES10 <= -1 )
replace VAGINDR10 = . if (VAGINDR10 <= -1 )
replace FEELBLU10 = . if (FEELBLU10 <= -1 )
replace DIZZY10 = . if (DIZZY10 <= -1 )
replace FORGET10 = . if (FORGET10 <= -1 )
replace MOODCHG10 = . if (MOODCHG10 <= -1 )
replace HARTRAC10 = . if (HARTRAC10 <= -1 )
replace FEARFULA10 = . if (FEARFULA10 <= -1 )
replace HDACHE10 = . if (HDACHE10 <= -1 )
replace BRSTPAI10 = . if (BRSTPAI10 <= -1 )
replace VAGIRRT10 = . if (VAGIRRT10 <= -1 )
replace VAGDISH10 = . if (VAGDISH10 <= -1 )
replace VAGSORE10 = . if (VAGSORE10 <= -1 )
replace TRBLSLE10 = . if (TRBLSLE10 <= -1 )
replace WAKEUP10 = . if (WAKEUP10 <= -1 )
replace WAKEARL10 = . if (WAKEARL10 <= -1 )
replace SLEEPQL10 = . if (SLEEPQL10 <= -1 )
replace INVOLEA10 = . if (INVOLEA10 <= -1 )
replace DAYSLEA10 = . if (DAYSLEA10 <= -1 )
replace COUGHLE10 = . if (COUGHLE10 <= -1 )
replace URGEVOI10 = . if (URGEVOI10 <= -1 )
replace OTHRLEA10 = . if (OTHRLEA10 <= -1 )
replace LEKAMNT10 = . if (LEKAMNT10 <= -1 )
replace LEKBOTH10 = . if (LEKBOTH10 <= -1 )
replace NUTRIRE10 = . if (NUTRIRE10 <= -1 )
replace HERBREM10 = . if (HERBREM10 <= -1 )
replace PSYCMET10 = . if (PSYCMET10 <= -1 )
replace PHYSMET10 = . if (PHYSMET10 <= -1 )
replace FOLKMED10 = . if (FOLKMED10 <= -1 )
replace OTHRTHE10 = . if (OTHRTHE10 <= -1 )
replace STARTNE10 = . if (STARTNE10 <= -1 )
replace WORKTRB10 = . if (WORKTRB10 <= -1 )
replace QUITJOB10 = . if (QUITJOB10 <= -1 )
replace WORKLOA10 = . if (WORKLOA10 <= -1 )
replace PRTUNEM10 = . if (PRTUNEM10 <= -1 )
replace MONEYPR10 = . if (MONEYPR10 <= -1 )
replace WORSREL10 = . if (WORSREL10 <= -1 )
replace RELATEN10 = . if (RELATEN10 <= -1 )
replace SERIPRO10 = . if (SERIPRO10 <= -1 )
replace CHILDMO10 = . if (CHILDMO10 <= -1 )
replace RESPCAR10 = . if (RESPCAR10 <= -1 )
replace LEGALPR10 = . if (LEGALPR10 <= -1 )
replace CRELDIE10 = . if (CRELDIE10 <= -1 )
replace CLOSDIE10 = . if (CLOSDIE10 <= -1 )
replace SELFVIO10 = . if (SELFVIO10 <= -1 )
replace FAMLVIO10 = . if (FAMLVIO10 <= -1 )
replace PHYSILL10 = . if (PHYSILL10 <= -1 )
replace MAJEVEN10 = . if (MAJEVEN10 <= -1 )
replace CRNTMAR10 = . if (CRNTMAR10 <= -1 )
replace RWRDREL10 = . if (RWRDREL10 <= -1 )
replace STRSREL10 = . if (STRSREL10 <= -1 )
replace INCOME10 = . if (INCOME10 <= -1 )
replace HOW_HAR10 = . if (HOW_HAR10 <= -1 )
replace INTRPAN10 = . if (INTRPAN10 <= -1 )
replace DISIPAN10 = . if (DISIPAN10 <= -1 )
replace EXCIPAN10 = . if (EXCIPAN10 <= -1 )
replace UPSEPAN10 = . if (UPSEPAN10 <= -1 )
replace STROPAN10 = . if (STROPAN10 <= -1 )
replace GUILPAN10 = . if (GUILPAN10 <= -1 )
replace SCARPAN10 = . if (SCARPAN10 <= -1 )
replace HOSTPAN10 = . if (HOSTPAN10 <= -1 )
replace ENTHPAN10 = . if (ENTHPAN10 <= -1 )
replace PROUPAN10 = . if (PROUPAN10 <= -1 )
replace IRRIPAN10 = . if (IRRIPAN10 <= -1 )
replace ALERPAN10 = . if (ALERPAN10 <= -1 )
replace ASHAPAN10 = . if (ASHAPAN10 <= -1 )
replace INSPPAN10 = . if (INSPPAN10 <= -1 )
replace NERVPAN10 = . if (NERVPAN10 <= -1 )
replace DETEPAN10 = . if (DETEPAN10 <= -1 )
replace ATTEPAN10 = . if (ATTEPAN10 <= -1 )
replace JITTPAN10 = . if (JITTPAN10 <= -1 )
replace ACTIPAN10 = . if (ACTIPAN10 <= -1 )
replace AFRAPAN10 = . if (AFRAPAN10 <= -1 )
replace COURTES10 = . if (COURTES10 <= -1 )
replace RESPECT10 = . if (RESPECT10 <= -1 )
replace POORSER10 = . if (POORSER10 <= -1 )
replace NOTSMAR10 = . if (NOTSMAR10 <= -1 )
replace AFRAIDO10 = . if (AFRAIDO10 <= -1 )
replace DISHONS10 = . if (DISHONS10 <= -1 )
replace BETTER10 = . if (BETTER10 <= -1 )
replace INSULTE10 = . if (INSULTE10 <= -1 )
replace HARASSE10 = . if (HARASSE10 <= -1 )
replace IGNORED10 = . if (IGNORED10 <= -1 )
replace BCRACE10 = . if (BCRACE10 <= -1 )
replace BCETHN10 = . if (BCETHN10 <= -1 )
replace BCGENDR10 = . if (BCGENDR10 <= -1 )
replace BCAGE10 = . if (BCAGE10 <= -1 )
replace BCINCML10 = . if (BCINCML10 <= -1 )
replace BCLANG10 = . if (BCLANG10 <= -1 )
replace BCWGHT10 = . if (BCWGHT10 <= -1 )
replace BCPHAPP10 = . if (BCPHAPP10 <= -1 )
replace BCORIEN10 = . if (BCORIEN10 <= -1 )
replace OTHEREX10 = . if (OTHEREX10 <= -1 )
replace FLGSAAV10 = . if (FLGSAAV10 <= -1 )
replace LANGSAB10 = . if (LANGSAB10 <= -1 )
replace IMPORSE10 = . if (IMPORSE10 <= -1 )
replace DESIRSE10 = . if (DESIRSE10 <= -1 )
replace ENGAGSE10 = . if (ENGAGSE10 <= -1 )
replace NOPARTN10 = . if (NOPARTN10 <= -1 )
replace PARTPRO10 = . if (PARTPRO10 <= -1 )
replace PHYSPRO10 = . if (PHYSPRO10 <= -1 )
replace NOSEXOT10 = . if (NOSEXOT10 <= -1 )
replace SATISFY10 = . if (SATISFY10 <= -1 )
replace KISSING10 = . if (KISSING10 <= -1 )
replace TOUCHIN10 = . if (TOUCHIN10 <= -1 )
replace ORALSEX10 = . if (ORALSEX10 <= -1 )
replace INTCOUR10 = . if (INTCOUR10 <= -1 )
replace AROUSED10 = . if (AROUSED10 <= -1 )
replace PELVIC10 = . if (PELVIC10 <= -1 )
replace LUBRICN10 = . if (LUBRICN10 <= -1 )
replace ABLECLM10 = . if (ABLECLM10 <= -1 )
replace IMPCLMX10 = . if (IMPCLMX10 <= -1 )
replace SATISFD10 = . if (SATISFD10 <= -1 )
replace FREQUEN10 = . if (FREQUEN10 <= -1 )
replace MEN6MOS10 = . if (MEN6MOS10 <= -1 )
replace MASTURB10 = . if (MASTURB10 <= -1 )
replace HTMETHO10 = . if (HTMETHO10 <= -1 )
replace HTSELF10 = . if (HTSELF10 <= -1 )
replace SCALE10 = . if (SCALE10 <= -1 )
replace WTSELF10 = . if (WTSELF10 <= -1 )
replace WASTMEA10 = . if (WASTMEA10 <= -1 )
replace HIPMEAS10 = . if (HIPMEAS10 <= -1 )
replace LANGCOG10 = . if (LANGCOG10 <= -1 )
replace IMEDTHR10 = . if (IMEDTHR10 <= -1 )
replace IMEDCH110 = . if (IMEDCH110 <= -1 )
replace IMEDHOU10 = . if (IMEDHOU10 <= -1 )
replace IMEDFIR10 = . if (IMEDFIR10 <= -1 )
replace IMEDFMN10 = . if (IMEDFMN10 <= -1 )
replace IMEDCLM10 = . if (IMEDCLM10 <= -1 )
replace IMEDCH210 = . if (IMEDCH210 <= -1 )
replace IMEDRES10 = . if (IMEDRES10 <= -1 )
replace IMEDMIN10 = . if (IMEDMIN10 <= -1 )
replace IMEDINJ10 = . if (IMEDINJ10 <= -1 )
replace IMEDEVR10 = . if (IMEDEVR10 <= -1 )
replace IMEDWEL10 = . if (IMEDWEL10 <= -1 )
replace SDMTSTA10 = . if (SDMTSTA10 <= -1 )
replace DIGIT1A10 = . if (DIGIT1A10 <= -1 )
replace DIGIT1B10 = . if (DIGIT1B10 <= -1 )
replace DIGIT2A10 = . if (DIGIT2A10 <= -1 )
replace DIGIT2B10 = . if (DIGIT2B10 <= -1 )
replace DIGIT3A10 = . if (DIGIT3A10 <= -1 )
replace DIGIT3B10 = . if (DIGIT3B10 <= -1 )
replace DIGIT4A10 = . if (DIGIT4A10 <= -1 )
replace DIGIT4B10 = . if (DIGIT4B10 <= -1 )
replace DIGIT5A10 = . if (DIGIT5A10 <= -1 )
replace DIGIT5B10 = . if (DIGIT5B10 <= -1 )
replace DIGIT6A10 = . if (DIGIT6A10 <= -1 )
replace DIGIT6B10 = . if (DIGIT6B10 <= -1 )
replace DLAYTHR10 = . if (DLAYTHR10 <= -1 )
replace DLAYCH110 = . if (DLAYCH110 <= -1 )
replace DLAYHOU10 = . if (DLAYHOU10 <= -1 )
replace DLAYFIR10 = . if (DLAYFIR10 <= -1 )
replace DLAYFMN10 = . if (DLAYFMN10 <= -1 )
replace DLAYCLM10 = . if (DLAYCLM10 <= -1 )
replace DLAYCH210 = . if (DLAYCH210 <= -1 )
replace DLAYRES10 = . if (DLAYRES10 <= -1 )
replace DLAYMIN10 = . if (DLAYMIN10 <= -1 )
replace DLAYINJ10 = . if (DLAYINJ10 <= -1 )
replace DLAYEVR10 = . if (DLAYEVR10 <= -1 )
replace DLAYWEL10 = . if (DLAYWEL10 <= -1 )
replace FLGCOG10 = . if (FLGCOG10 <= -1 )
replace BIODAY10 = . if (BIODAY10 <= -1 )
replace LANGBIO10 = . if (LANGBIO10 <= -1 )
replace COMPBIA10 = . if (COMPBIA10 <= -1 )
replace AICDPUM10 = . if (AICDPUM10 <= -1 )
replace EXER12H10 = . if (EXER12H10 <= -1 )
replace EAT5HR10 = . if (EAT5HR10 <= -1 )
replace ALCO24H10 = . if (ALCO24H10 <= -1 )
replace EMBDDEV10 = . if (EMBDDEV10 <= -1 )
replace METJEWL10 = . if (METJEWL10 <= -1 )
replace ONMEASS10 = . if (ONMEASS10 <= -1 )
replace SIDE10 = . if (SIDE10 <= -1 )
replace CONDRAW10 = . if (CONDRAW10 <= -1 )
replace IMPERAW10 = . if (IMPERAW10 <= -1 )
replace BIORRUN10 = . if (BIORRUN10 <= -1 )
replace BLDDRAW10 = . if (BLDDRAW10 <= -1 )
replace PREGNAN10 = . if (PREGNAN10 <= -1 )
replace EATDRIN10 = . if (EATDRIN10 <= -1 )
replace STRTPER10 = . if (STRTPER10 <= -1 )
replace BLDRWAT10 = . if (BLDRWAT10 <= -1 )
replace ALCHL2410 = . if (ALCHL2410 <= -1 )
replace FLGCV10 = . if (FLGCV10 <= -1 )
replace FLGDIF10 = . if (FLGDIF10 <= -1 )
replace FLGHRMV10 = . if (FLGHRMV10 <= -1 )
replace SPSCTIM10 = "" if (SPSCTIM10 == ".")
replace HPSCTIM10 = "" if (HPSCTIM10 == ".")
replace BMDFLG10 = . if (BMDFLG10 <= -1 )
replace MISSCON10 = . if (MISSCON10 <= -1 )
replace MISSPHY10 = . if (MISSPHY10 <= -1 )
replace FLAGSRP10 = . if (FLAGSRP10 <= -1 )
replace FLGBIOV10 = . if (FLGBIOV10 <= -1 )


