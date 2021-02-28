#!/bin/bash
prody catdcd ./3_WT/3_production/seokWT_100to250ns.dcd --psf ./3_WT/3_production/seokWT.psf -o ./3_WT/3_production/seokWT_100to250ns_each25k.dcd -s protein --stride 5
prody catdcd ./3_WT/4_production_repeat2/seokWT_100to250ns_r2.dcd --psf ./3_WT/3_production/seokWT.psf -o ./3_WT/4_production_repeat2/seokWT_100to250ns_r2_each25k.dcd -s protein --stride 5

prody catdcd ./4_P389L/3_production/seokP389L_100to250ns.dcd --psf ./4_P389L/3_production/seokP389L.psf -o ./4_P389L/3_production/seokP389L_100to250ns_each25k.dcd -s protein --stride 5
prody catdcd ./4_P389L/4_production_repeat2/seokP389L_100to250ns_r2.dcd --psf ./4_P389L/3_production/seokP389L.psf -o ./4_P389L/4_production_repeat2/seokP389L_100to250ns_r2_each25k.dcd -s protein --stride 5

prody catdcd ./5_G317R/3_production/seokG317R_100to250ns.dcd --psf ./5_G317R/3_production/seokG317R.psf -o ./5_G317R/3_production/seokG317R_100to250ns_each25k.dcd -s protein --stride 5
prody catdcd ./5_G317R/4_production_repeat2/seokG317R_100to250ns_r2.dcd --psf ./5_G317R/3_production/seokG317R.psf -o ./5_G317R/4_production_repeat2/seokG317R_100to250ns_r2_each25k.dcd -s protein --stride 5

prody catdcd ./6_L430P/3_production/seokL430P_100to250ns.dcd --psf ./6_L430P/3_production/seokL430P.psf -o ./6_L430P/3_production/seokL430P_100to250ns_each25k.dcd -s protein --stride 5
prody catdcd ./6_L430P/4_production_repeat2/seokL430P_100to250ns_r2.dcd --psf ./6_L430P/3_production/seokL430P.psf -o ./6_L430P/4_production_repeat2/seokL430P_100to250ns_r2_each25k.dcd -s protein --stride 5

prody catdcd ./7_S50P/3_production/seokS50P_100to250ns.dcd --psf ./7_S50P/3_production/seokS50P.psf -o ./7_S50P/3_production/seokS50P_100to250ns_each25k.dcd -s protein --stride 5
prody catdcd ./7_S50P/4_production_repeat2/seokS50P_100to250ns_r2.dcd --psf ./7_S50P/3_production/seokS50P.psf -o ./7_S50P/4_production_repeat2/seokS50P_100to250ns_r2_each25k.dcd -s protein --stride 5

prody catdcd ./8_R266Q/3_production/seokR266Q_100to250ns.dcd --psf ./8_R266Q/3_production/seokR266Q.psf -o ./8_R266Q/3_production/seokR266Q_100to250ns_each25k.dcd -s protein --stride 5
prody catdcd ./8_R266Q/4_production_repeat2/seokR266Q_100to250ns_r2.dcd --psf ./8_R266Q/3_production/seokR266Q.psf -o ./8_R266Q/4_production_repeat2/seokR266Q_100to250ns_r2_each25k.dcd -s protein --stride 5
