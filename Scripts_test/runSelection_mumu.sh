# Pre
./loadRoot.sh

# Jobs
./runHandler.sh     mumu    muon_2016       1
./runHandler.sh     mumu    muon_2016       2
./runHandler.sh     mumu    muon_2016       3
./runHandler.sh     mumu    muon_2016       4
./runHandler.sh     mumu    muon_2016       5
./runHadd.sh        mumu    muon_2016

./runHandler.sh     mumu    zjets_m-10to50  1
./runHadd.sh        mumu    zjets_m-10to50

./runHandler.sh     mumu    zjets_m-50      1
./runHandler.sh     mumu    zjets_m-50      2
./runHadd.sh        mumu    zjets_m-50 

#/runHandler.sh     mumu    dy_m-10to50     1
#/runHadd.sh        mumu    dy_m-10to50

#/runHandler.sh     mumu    dy_m-50         1
#/runHandler.sh     mumu    dy_m-50         2
#/runHandler.sh     mumu    dy_m-50         3
#/runHandler.sh     mumu    dy_m-50         4
#/runHandler.sh     mumu    dy_m-50         5
#/runHadd.sh        mumu    dy_m-50

./runHandler.sh     mumu    ttbar           1
./runHadd.sh        mumu    ttbar

./runHandler.sh     mumu    ttz_2l2nu       1
./runHadd.sh        mumu    ttz_2l2nu

./runHandler.sh     mumu    ww_2l2nu        1
./runHadd.sh        mumu    ww_2l2nu

./runHandler.sh     mumu    wz_3lnu         1
./runHadd.sh        mumu    wz_3lnu

./runHandler.sh     mumu    zz_4l           1
./runHadd.sh        mumu    zz_4l

# Post
./runStacker.sh     mumu
