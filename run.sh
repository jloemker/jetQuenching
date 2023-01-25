export OPTIONS="-b --configuration json://config.json"

o2-analysistutorial-correlationv0jet ${OPTIONS} |  o2-analysis-timestamp ${OPTIONS} | o2-analysis-event-selection ${OPTIONS} | o2-analysis-trackselection ${OPTIONS} | o2-analysis-trackextension ${OPTIONS} | o2-analysis-lf-lambdakzerobuilder ${OPTIONS} | o2-analysis-pid-tpc ${OPTIONS} | o2-analysis-pid-tof-base ${OPTIONS} | o2-analysis-pid-tof ${OPTIONS} | o2-analysis-multiplicity-table ${OPTIONS} 


#o2-analysistutorial-correlationv0jet ${OPTIONS} |  o2-analysis-timestamp ${OPTIONS} | o2-analysis-event-selection ${OPTIONS} | o2-analysis-trackselection ${OPTIONS} | o2-analysis-trackextension ${OPTIONS} | o2-analysis-lf-lambdakzerobuilder ${OPTIONS} | o2-analysis-pid-tpc ${OPTIONS} | o2-analysis-pid-tof-base ${OPTIONS} | o2-analysis-pid-tof ${OPTIONS} | o2-analysis-multiplicity-table ${OPTIONS} | o2-analysis-mc-converter ${OPTIONS}
 #| o2-analysis-lf-lambdakzerobuilder ${OPTIONS}

#o2-analysistutorial-correlationv0jet ${OPTIONS} |  o2-analysis-timestamp ${OPTIONS} | o2-analysis-event-selection ${OPTIONS} | o2-analysis-je-filter ${OPTIONS} | o2-analysis-trackselection ${OPTIONS} | o2-analysis-trackextension ${OPTIONS} | o2-analysis-lf-lambdakzerobuilder ${OPTIONS} | o2-analysis-pid-tpc ${OPTIONS} | o2-analysis-multiplicity-table ${OPTIONS}

#o2-analysis-je-filter ${OPTIONS} | o2-analysis-timestamp ${OPTIONS} | o2-analysis-event-selection ${OPTIONS} | o2-analysis-trackselection ${OPTIONS} | o2-analysis-trackextension ${OPTIONS}| o2-analysis-lf-lambdakzerobuilder ${OPTIONS} | o2-analysis-pid-tpc ${OPTIONS} | o2-analysis-multiplicity-table ${OPTIONS} | o2-analysistutorial-correlationv0jet ${OPTIONS} 
#| o2-analysis-lf-lambdakzerobuilder ${OPTIONS} gives me [ERROR] SEVERE: Device lambdakzero-builder (197037) had at least one message above severity 4: Unhandled o2::framework::runtime_error reached the top of main of o2-analysis-lf-lambdakzerobuilder, device shutting down. Reason: invalid track covariance
#| o2-analysis-lf-lambdakzerofinder ${OPTIONS}
#| o2-analysis-track-propagation ${OPTIONS}
