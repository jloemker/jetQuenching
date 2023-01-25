export OPTIONS="-b --configuration json://configAO2D.json"

o2-analysistutorial-correlationv0jet ${OPTIONS} |  o2-analysis-timestamp ${OPTIONS} | o2-analysis-event-selection ${OPTIONS} | o2-analysis-trackselection ${OPTIONS} | o2-analysis-trackextension ${OPTIONS} | o2-analysis-lf-lambdakzerobuilder ${OPTIONS} | o2-analysis-pid-tpc ${OPTIONS} | o2-analysis-pid-tof-base ${OPTIONS} | o2-analysis-pid-tof ${OPTIONS} | o2-analysis-multiplicity-table ${OPTIONS} 

