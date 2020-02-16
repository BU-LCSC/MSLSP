



planted_2018 = min(progress$Week.Ending[progress$State = 'MISSISSIPPI' & progress$Year == 2018 & progress$Data.Item == "CORN - PROGRESS, MEASURED IN PCT PLANTED" & progress$Value > 50])
emerged_2018 = min(progress$Week.Ending[progress$State = 'MISSISSIPPI' & progress$Year == 2018 & progress$Data.Item == "CORN - PROGRESS, MEASURED IN PCT EMERGED" & progress$Value > 50])
silking_2018 = min(progress$Week.Ending[progress$State = 'MISSISSIPPI' & progress$Year == 2018 & progress$Data.Item == "CORN - PROGRESS, MEASURED IN PCT SILKING" & progress$Value > 50])
dough_2018 = min(progress$Week.Ending[progress$State = 'MISSISSIPPI' & progress$Year == 2018 & progress$Data.Item == "CORN - PROGRESS, MEASURED IN PCT DOUGH" & progress$Value > 50])
dent_2018 = min(progress$Week.Ending[progress$State = 'MISSISSIPPI' & progress$Year == 2018 & progress$Data.Item == "CORN - PROGRESS, MEASURED IN PCT DENTED" & progress$Value > 50])
mature_2018 = min(progress$Week.Ending[progress$State = 'MISSISSIPPI' & progress$Year == 2018 & progress$Data.Item == "CORN - PROGRESS, MEASURED IN PCT MATURE" & progress$Value > 50])
harvested_2018 = min(progress$Week.Ending[progress$State = 'MISSISSIPPI' & progress$Year == 2018 & progress$Data.Item == "CORN, GRAIN - PROGRESS, MEASURED IN PCT HARVESTED" & progress$Value > 50])