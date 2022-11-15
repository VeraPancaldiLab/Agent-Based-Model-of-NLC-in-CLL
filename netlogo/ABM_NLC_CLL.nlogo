;;;;;;;;;;;;;;;;;;
;;;;; BREEDS ;;;;;
;;;;;;;;;;;;;;;;;;

breed [ NeedSignalCancerCells NeedSignalCancerCell ]
breed [ ApoCancerCells ApoCancerCell ]
breed [ DeadCancerCells DeadCancerCell ]
breed [ Monocytes Monocyte ]
breed [ NurseLikeCells NurseLikeCell ]
breed [ Macrophages Macrophage ]


;;;;;;;;;;;;;;;;;;;;;
;;;;; VARIABLES ;;;;;
;;;;;;;;;;;;;;;;;;;;;

globals [

  ;; create limits of the world
  ;z-max ; for the 3D version of the model
  coef

  ;; global variables for headless mode
  simulation-duration
  density
  nb-cancer-cells-init ; initial number of cancer cells
  prop-monocytes-init ; initial proportion of monocytes
  prop-apo-init ; initial proportion of apoptotic cancer cells
  apoptosis-threshold ; threshold from need-signal to apoptotic state
  death-threshold ; threshold from apoptotic to dead state
  my-seed

  ;; Parameters to optimize
  gamma-life-init-shape ; Shape parameter of the Gamma distribution setting NeedSignal cell life at initialization (alpha)
  gamma-life-init-rate ; Rate parameter of the Gamma distribution setting NeedSignal cell life at initialization (beta)
  need-sig-cells-mvt-proba ; NeedSignal cell movement probability
  cll-sensing-distance ; NeedSignal cell sensing distance to NLCs
  apo-cells-movement-proba ; Apoptotic cell movement probability
  mono-phago-eff ; Monocyte phagocytosis efficiency
  adherence-init-std ; Standard deviation parameter of the Normal distribution setting Monocyte differentiation status at initialization
  adherence-time ; Monocyte differentiation threshold
  monocyte-sensing-distance ; Monocyte sensing distance to Dead and Apoptotic cancer cells
  macro-phago-eff ; Macrophage phagocytosis efficiency
  macro-kill-eff ; Macrophage killing efficiency
  signal-init-mean ; Mean parameter of the Normal distribution setting Macrophage polarization status at instantiation
  signal-init-std ; Standard deviation parameter of the Normal distribution setting Macrophage polarization status at instantiation
  macrophage-sensing-distance ; Macrophage sensing distance to Dead and Apoptotic cancer cells
  nlc-phago-eff ; NLC phagocytosis efficiency
  anti-apo-boost ; Level of protective effect from the anti apoptotic signals
  layers-around-NLC ; Number of layers around NLCs
  nlc-sensing-distance ; NLC sensing distance to Dead and Apoptotic cancer cells
  nlc-threshold ; NLC threshold

  ;; outputs for headless mode
  nb-NeedSignalCancerCells
  nb-ApoCancerCells
  nb-Monocytes
  nb-NurseLikeCells
  nb-DeadCancerCells
  nb-monocytes-init

  ;; lists for time series in headless mode
  ts-Viability
  ts-Concentration
]

patches-own [
  forbidden ;; manual fix to prevent auto-wrapping
  Chemical_A ;; anti-apoptotic chemokine secreted by NurseLikeCells
]

NeedSignalCancerCells-own[
  Life ;; the lower the Life, the more apoptotic the cell is. When reaching 0, the cell enters apoptosis which is irreversible.
]

ApoCancerCells-own[
  Life
]

DeadCancerCells-own[
  Life
]

Monocytes-own[
  NLC-polarization ;; NLC-polarization signal. It can be seen as a proxy for the polarization level of monocytes into NLCs.
  Adherence ;; time spent after cell seeding
  counter-contacts
  contacts
]

NurseLikeCells-own[
  NLC-polarization
  counter-contacts
  contacts
]

Macrophages-own[
  NLC-polarization
  counter-contacts
  contacts
]

;;;;;;;;;;;;;;;;;
;;;;; SETUP ;;;;;
;;;;;;;;;;;;;;;;;

to setup
  clear-all
  reset-ticks
  set my-seed random 1000000
  setup-Globals
  setup-World
  setup-Monocytes
  setup-NeedSignalCancerCells
  setup-ApoCancerCells
end

to headless-setup
  reset-ticks
  set my-seed random 1000000
  setup-World
  setup-Monocytes
  setup-NeedSignalCancerCells
  setup-ApoCancerCells
  setup-ts
end

to setup-Globals
; fixed
  ;set z-max 0
  set coef 1
  set apoptosis-threshold 0
  set death-threshold -500
  set density 55
  set nb-cancer-cells-init gui-nb-cancer-cells-init
  set prop-monocytes-init gui-prop-mono-init
  set prop-apo-init gui-prop-apo-init

; optimized
  set gamma-life-init-rate gui-gamma-life-init-rate
  set gamma-life-init-shape gui-life-init-shape
  set nlc-threshold gui-nlc-threshold
  set anti-apo-boost gui-anti-apo-boost
  set mono-phago-eff gui-mono-phago-eff
  set nlc-phago-eff gui-nlc-phago-eff
  set macro-phago-eff gui-macro-phago-eff
  set macro-kill-eff gui-macro-kill-eff
  set adherence-time gui-diff-mean
  set adherence-init-std gui-diff-std
  set signal-init-mean gui-sig-init
  set signal-init-std gui-sig-init-std
  set apo-cells-movement-proba gui-apo-mov
  set need-sig-cells-mvt-proba gui-need-sig-mov
  set layers-around-NLC gui-layers
  set cll-sensing-distance gui-cll-sens-dist
  set monocyte-sensing-distance gui-mono-sens-dist
  set nlc-sensing-distance gui-nlc-sens-dist
  set macrophage-sensing-distance gui-macro-sens-dist

end

to setup-World
  ;; get initial number of cells
  set nb-monocytes-init ceiling (prop-monocytes-init * nb-cancer-cells-init / 100)
  let nb-cells-init nb-cancer-cells-init + nb-monocytes-init

;;; use these lines to make the model work in 3D
  ;let world-radius ceiling sqrt (nb-cells-init * 100 / (density * pi * (z-max * 2 + 1)))
  ;resize-world (- world-radius * coef) (world-radius * coef) (- world-radius * coef) (world-radius * coef) (- z-max * coef) (z-max * coef)

  let world-radius ceiling sqrt (nb-cells-init * 100 / (density * pi ))
  resize-world (- world-radius * coef) (world-radius * coef) (- world-radius * coef) (world-radius * coef)

  ;; close the world
  ;; only not forbidden patches are accessible to the agents
  ask patches [
    ;ifelse ((pxcor * pxcor + pycor * pycor) > (world-radius * world-radius) or abs pzcor > z-max) ; for the 3D version
    ifelse ((pxcor * pxcor + pycor * pycor) > (world-radius * world-radius) )
    [ set forbidden true ]
    [ set forbidden false
      set Chemical_A 0
    ]
  ]

   ;; create NeedSignalCancerCells
  ask n-of nb-cells-init (patches with [not forbidden]) [ sprout-NeedSignalCancerCells 1 ]


end

to setup-Monocytes
  ;; create Monocytes
  let myMonocytes n-of nb-monocytes-init turtles
  ask myMonocytes
  [
    set breed Monocytes
    set NLC-polarization random-normal signal-init-mean signal-init-std
    set Adherence random-normal (0) (adherence-init-std)
    set counter-contacts 0
    set contacts 0
    set color blue
    set shape "pentagon"
    set size 1
  ]
end

to setup-NeedSignalCancerCells
  ask NeedSignalCancerCells
  [
    let my_life_init_shape gamma-life-init-shape
    let my_life_init_scale 1 / gamma-life-init-rate
    set Life random-gamma my_life_init_shape my_life_init_scale
    set color red
    if Life <= 0
    [
      set breed ApoCancerCells
      set color yellow
    ]
    set shape "default"
    set size 1
  ]
end

to setup-ApoCancerCells
  ;; create ApoCancerCells
  let myApoCancerCells n-of ( nb-cancer-cells-init * prop-apo-init / 100 ) turtles
  ask myApoCancerCells
  [
    set breed ApoCancerCells
    set Life -1
    set color yellow
    set shape "default"
    set size 1
  ]
end

to setup-ts
  ;; update time series
  set ts-Viability (list (getViability))
  set ts-Concentration (list (getRemainingCellRatio))
end

;;;;;;;;;;;;;;
;;;;; UTILS ;;
;;;;;;;;;;;;;;

to-report getSeed
  report my-seed
end

to-report getNumberNLC
  report count NurseLikeCells
end

to-report getNLCFraction
  report (count NurseLikeCells / nb-monocytes-init) * 100
end

to-report getMacroFraction
  report (count Macrophages / nb-monocytes-init) * 100
end

to-report getMonoFraction
  report (count Monocytes / nb-monocytes-init) * 100
end

to-report getNumberMono
  report count Monocytes
end

to-report getNumberMacrophages
  report count Macrophages
end

to-report getDead
  report count DeadCancerCells
end

to-report getNeedSignal
  report count NeedSignalCancerCells
end

to-report getLateApo
  report count ApoCancerCells
end

to-report getdeath-threshold
  report death-threshold
end

to-report getViability ;; cancer B-CLL cell viability
  report ((count NeedSignalCancerCells ) / (count DeadCancerCells + count ApoCancerCells + count NeedSignalCancerCells )) * 100
end

to-report getRemainingCellRatio
  report ((count DeadCancerCells + count ApoCancerCells + count NeedSignalCancerCells) / (nb-cancer-cells-init)) * 100
end

to-report getCounterContacts
  report (counter-contacts)
end

to-report getContacts
  report contacts
end

;;;;;;;;;;;;;;;;;;;;;
;;;;; ACTIONS ;;;;;;;
;;;;;;;;;;;;;;;;;;;;;

to do_phagocytosis [ phago-eff sensing-distance ]
  ;; look for dead or apoptotic cells in the neighbors
  let debris (turtle-set (DeadCancerCells-on neighbors) (ApoCancerCells-on neighbors))
  ifelse (any? debris)
  [ if (random 100 < phago-eff) [ ask one-of debris [die] ] ]
  [ move sensing-distance nobody ]
end

to kill_CLL [ kill-eff sensing-distance ]
  ;; look for need signal CLL cells in the neighbors
  let CLL (turtle-set (NeedSignalCancerCells-on neighbors))
  ifelse (any? CLL)
  [if (random 100 < kill-eff) [ ask one-of CLL [die] ] ]
  [ move sensing-distance nobody ]
end

to move [steps goal] ;; agent procedure
  ifelse (goal != nobody)
  ;; move towards goal
  [
    while [steps > 0] [
      ;; distance from the agent to the goal
      let current-distance distance goal
      ;; find neighbors who are closer to the goal
      let closer-positions neighbors with [not forbidden and
                                           not any? turtles-here and
                                           distance goal <= current-distance and
                                           distance goal > 1]
      ifelse (any? closer-positions)
      ;; move closer to the goal
      [
        move-to min-one-of closer-positions [distance goal]
        set steps steps - 1
      ]
      ;; impossible to move closer, stay put
      [ set steps 0 ]
    ]
  ]
  ;; random move
  [
    while [steps > 0] [
      ;; find accessible positions
      let possible-positions neighbors with [not forbidden and not any? turtles-here]
      ifelse (any? possible-positions)
      ;; move
      [
        move-to one-of possible-positions
        set steps steps - 1
      ]
      ;; impossible to move closer, stay put
      [ set steps 0 ]
    ]
  ]
end


to move-to-apoptotic [ my_radius ] ;; myeloid cell procedure
  let possible-positions patches in-radius my_radius with [not forbidden and not any? turtles-here and any? ApoCancerCells-on neighbors or any? DeadCancerCells-on neighbors]
  if (any? possible-positions)
  [ move-to one-of possible-positions ]
end

to move-to-NLC [sensing-distance layers] ;; cancer cell procedure
  ;; go where an NLC is adjacent (or distant of n=layers patches). If there is no NLC move randomly (or stay put)
  let possible-positions patches in-radius sensing-distance with [not forbidden and not any? turtles-here and any? NurseLikeCells in-radius layers]
  ifelse (any? possible-positions)
  [ if (random 100 < need-sig-cells-mvt-proba * 10)[ move-to one-of possible-positions ] ]
  [ if (random 100 < need-sig-cells-mvt-proba * 10)[ move sensing-distance nobody ] ]
end

;;;;;;;;;;;;;;
;;;;; GO ;;;;;
;;;;;;;;;;;;;;

to go
  ;reset-timer
  ifelse (count NeedSignalCancerCells + count ApoCancerCells  > 0)
  [
    common-go-procedures
    tick
  ]
  [ stop ]
end

to headless-go
  common-go-procedures
  ;; update time series of cell counts
  save-ts
  ;; increment the time step
  tick
end

to common-go-procedures ;; between headless and gui modes
  update-positions
  update-phagocytosis
  update-chemicals
  update-breeds
end

to update-positions

  ;; NeedSignalCancerCells actions
  ask NeedSignalCancerCells
  [ move-to-NLC cll-sensing-distance layers-around-NLC ]

  ;; ApoCancerCells actions
  ask ApoCancerCells
  [ if ( random 100 < apo-cells-movement-proba * 10 ) [move cll-sensing-distance nobody ]]

  ; Monocytes actions
  ask Monocytes
  [ move-to-apoptotic monocyte-sensing-distance ]

  ;; NurseLikeCells actions
  ask NurseLikeCells
  [ move-to-apoptotic nlc-sensing-distance ]

  ;; Macrophage actions
  ask Macrophages
  [ move-to-apoptotic macrophage-sensing-distance ]

end


to update-phagocytosis
  ; Monocytes actions
  ask Monocytes

  [
    do_phagocytosis mono-phago-eff monocyte-sensing-distance
  ]

  ;; NurseLikeCells actions
  ask NurseLikeCells
  [
    do_phagocytosis nlc-phago-eff nlc-sensing-distance
  ]

  ;; Macrophage actions
  ask Macrophages
  [
    kill_CLL macro-kill-eff macrophage-sensing-distance
    do_phagocytosis macro-phago-eff macrophage-sensing-distance
  ]
end


to update-chemicals

  ;; Patches actions
  diffuse chemical_A (1 / 1000)

  ;; NeedSignalCancerCells actions
  ask NeedSignalCancerCells
  [
    ifelse ( [Chemical_A] of patch-here >= 1)
    [
      set Life Life + anti-apo-boost;
      set Chemical_A Chemical_A - 1
    ]
    [ set Life Life - 1 ]
  ]

  ;; ApoCancerCells actions
  ask ApoCancerCells
  [
    set Life Life - 1
  ]

  ;; DeadCancerCells actions
  ask DeadCancerCells
  [
    set Life Life - 1
  ]

  ;; Monocytes actions
  ask Monocytes
  [
    set Adherence Adherence + 1
  ]

  ;; NurseLikeCells actions
  ask NurseLikeCells
  [
    ifelse (any? NeedSignalCancerCells in-radius 1.5)
    [ set Chemical_A Chemical_A + 1 ]
    [ set NLC-polarization NLC-polarization - 1 ]
  ]

  ;; Macrophages actions
  ask Macrophages
  [
    ifelse (any? NeedSignalCancerCells in-radius 1.5 or any? ApoCancerCells in-radius 1.5)
    [ set NLC-polarization NLC-polarization + 1 ]
    [ set NLC-polarization NLC-polarization - 1 ]
  ]


end

to update-breeds

  ;; NeedSignalCancerCells
  ask NeedSignalCancerCells
  [
    if ( Life <= apoptosis-threshold )
    [
      ;; change from NeedSignalCancerCell to ApoCancerCell
      set breed ApoCancerCells
      set color yellow
      set shape "default"
    ]
  ]

  ;; ApoCancerCells
  ask ApoCancerCells
  [
    if ( Life <= death-threshold )
    [
      ;; change from ApoCancerCell to DeadCancerCell
      set breed DeadCancerCells
      set color grey
      set shape "default"
    ]
  ]

  ; Monocytes
  ask Monocytes
  [
    if (Adherence >= adherence-time)
    [
      ;; change from Monocyte to Macrophage
        set breed Macrophages
        set color orange
        set shape "pentagon"
        set size 2
    ]
  ]

  ;; Macrophages
  ask Macrophages
  [
    if ( NLC-polarization >= nlc-threshold )
    [
       ;; change from Macropahges to NurseLikeCells
       set breed NurseLikeCells
       set color green
       set shape "pentagon"
       set size 2
    ]
  ]


end

to save-ts
  if (member? ticks [0 24 48 72 144 168 192 216 240 312]) ; saving time points to calculate fitness functions
  [
    set ts-Viability lput (getViability) ts-Viability
    set ts-Concentration lput (getRemainingCellRatio) ts-Concentration
  ]
end

to outputs
  ;; Outputs
  set nb-NeedSignalCancerCells count NeedSignalCancerCells
  set nb-ApoCancerCells count ApoCancerCells
  set nb-Monocytes count Monocytes
  set nb-NurseLikeCells  count NurseLikeCells
  set nb-DeadCancerCells count DeadCancerCells
  set nb-monocytes-init ceiling (prop-monocytes-init * nb-cancer-cells-init / 100)
end
@#$#@#$#@
GRAPHICS-WINDOW
750
10
1565
826
-1
-1
7.02
1
10
1
1
1
0
1
1
1
-57
57
-57
57
0
0
1
ticks
30.0

BUTTON
0
10
73
43
NIL
setup
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
74
43
161
76
go once
go
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

MONITOR
73
669
143
714
cancer cells
count ApoCancerCells + count NeedSignalCancerCells + count DeadCancerCells
17
1
11

TEXTBOX
217
13
334
34
model
15
15.0
1

MONITOR
140
669
203
714
monocytes
count Monocytes
17
1
11

MONITOR
278
669
328
714
NLCs
count NurseLikeCells
17
1
11

BUTTON
0
43
73
76
NIL
go
T
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

MONITOR
0
669
73
714
NIL
count turtles
17
1
11

MONITOR
155
758
261
803
allowed patches
count patches with [not forbidden]
17
1
11

PLOT
0
387
299
523
Cell populations
time
Cell types
0.0
10.0
0.0
10.0
true
true
"" ""
PENS
"dead" 1.0 0 -7500403 true "" "plot count DeadCancerCells"
"need-signal" 1.0 0 -2674135 true "" "plot count NeedSignalCancerCells"
"apoptotic" 1.0 0 -7171555 true "" "plot count ApoCancerCells"
"apo+dead" 1.0 0 -6459832 true "" "plot (count ApoCancerCells + count DeadCancerCells)"

MONITOR
0
713
73
758
dead cells
count DeadCancerCells
17
1
11

PLOT
0
521
300
666
monocytes
NIL
NIL
0.0
10.0
0.0
10.0
true
true
"" ""
PENS
"Monocytes" 1.0 0 -13345367 true "" "plot count Monocytes"
"NLCs" 1.0 0 -10899396 true "" "plot count NurseLikeCells"
"Macrophages" 1.0 0 -955883 true "" "plot count Macrophages "

MONITOR
138
713
205
758
apoptotic
count ApoCancerCells
17
1
11

MONITOR
399
758
475
803
days
ticks / 24
1
1
11

PLOT
299
521
548
666
Concentration
time
%
0.0
10.0
0.0
120.0
true
true
"" ""
PENS
"Concentration" 1.0 0 -2674135 true "" "if (nb-cancer-cells-init > 0)\n[plot ( ( (count DeadCancerCells + count ApoCancerCells + count NeedSignalCancerCells ) / nb-cancer-cells-init ) * 100)]"

MONITOR
411
669
475
714
Avg Signal
(sum ([NLC-polarization] of Monocytes) +  sum ([NLC-polarization] of Macrophages) + sum ([NLC-polarization] of NurseLikeCells))/ (count Monocytes + count Macrophages + count NurseLikeCells)
2
1
11

MONITOR
0
758
155
803
Chemokine on patches (avg)
(sum ([Chemical_A] of patches with [not forbidden]) )/ (count patches with [not forbidden])
2
1
11

PLOT
298
387
549
522
viability
time
viability
0.0
10.0
80.0
100.0
true
true
"" ""
PENS
"viability" 1.0 0 -16777216 true "" "if (count DeadCancerCells + count ApoCancerCells + count NeedSignalCancerCells) > 0\n[plot ( count NeedSignalCancerCells / (count DeadCancerCells + count ApoCancerCells + count NeedSignalCancerCells)) * 100]\n"

MONITOR
73
713
139
758
need-signal
count NeedSignalCancerCells
17
1
11

MONITOR
205
714
259
759
dead
count DeadCancerCells
17
1
11

MONITOR
202
669
279
714
Macrophages
count Macrophages
17
1
11

MONITOR
259
714
385
759
getCounterContacts
(sum ([counter-contacts] of monocytes) + sum ([counter-contacts] of Macrophages) )/ (count monocytes + count Macrophages) / ticks
17
1
11

MONITOR
385
714
475
759
nlc-threshold
nlc-threshold
17
1
11

MONITOR
261
758
333
803
Max Signal
max [NLC-polarization] of (turtle-set (Monocytes) (Macrophages) (NurseLikeCells))
17
1
11

MONITOR
332
758
399
803
min Signal
min [NLC-polarization] of (turtle-set (Monocytes) (Macrophages) (NurseLikeCells))
17
1
11

MONITOR
328
669
411
714
NLC Fraction
(count NurseLikeCells / nb-monocytes-init) * 100
2
1
11

INPUTBOX
0
131
101
191
gui-prop-mono-init
1.28
1
0
Number

INPUTBOX
450
267
561
327
gui-macro-phago-eff
92.0
1
0
Number

INPUTBOX
340
267
450
327
gui-mono-phago-eff
25.0
1
0
Number

INPUTBOX
560
267
656
327
gui-nlc-phago-eff
4.0
1
0
Number

INPUTBOX
487
327
569
387
gui-macro-kill-eff
0.0
1
0
Number

INPUTBOX
573
149
674
209
gui-anti-apo-boost
244.0
1
0
Number

INPUTBOX
283
191
353
251
gui-apo-mov
8.0
1
0
Number

INPUTBOX
191
191
284
251
gui-need-sig-mov
9.0
1
0
Number

INPUTBOX
424
327
488
387
gui-diff-std
6.0
1
0
Number

INPUTBOX
345
327
424
387
gui-diff-mean
49.0
1
0
Number

INPUTBOX
674
149
734
209
gui-layers
2.0
1
0
Number

INPUTBOX
354
208
439
268
gui-cll-sens-dist
3.0
1
0
Number

INPUTBOX
438
208
539
268
gui-mono-sens-dist
2.0
1
0
Number

INPUTBOX
644
208
734
268
gui-nlc-sens-dist
1.0
1
0
Number

INPUTBOX
539
208
644
268
gui-macro-sens-dist
2.0
1
0
Number

INPUTBOX
568
327
659
387
gui-nlc-threshold
99.0
1
0
Number

INPUTBOX
655
267
734
327
gui-sig-init
63.0
1
0
Number

INPUTBOX
657
327
734
387
gui-sig-init-std
4.0
1
0
Number

INPUTBOX
100
131
191
191
gui-prop-apo-init
4.55
1
0
Number

MONITOR
547
621
659
666
Concentration (%)
((count DeadCancerCells + count ApoCancerCells + count NeedSignalCancerCells ) / nb-cancer-cells-init ) * 100
2
1
11

MONITOR
547
577
659
622
Viability (%)
( count NeedSignalCancerCells / (count DeadCancerCells + count ApoCancerCells + count NeedSignalCancerCells)) * 100
2
1
11

INPUTBOX
190
72
312
132
gui-nb-cancer-cells-init
5000.0
1
0
Number

INPUTBOX
190
131
312
191
gui-gamma-life-init-rate
2120.0
1
0
Number

INPUTBOX
311
131
407
191
gui-life-init-shape
0.3188403442070939
1
0
Number

@#$#@#$#@
## WHAT IS IT?

(a general understanding of what the model is trying to show or explain)

## HOW IT WORKS

(what rules the agents use to create the overall behavior of the model)

## HOW TO USE IT

(how to use the model, including a description of each of the items in the Interface tab)

## THINGS TO NOTICE

(suggested things for the user to notice while running the model)

## THINGS TO TRY

(suggested things for the user to try to do (move sliders, switches, etc.) with the model)

## EXTENDING THE MODEL

(suggested things to add or change in the Code tab to make the model more complicated, detailed, accurate, etc.)

## NETLOGO FEATURES

(interesting or unusual features of NetLogo that the model uses, particularly in the Code tab; or where workarounds were needed for missing features)

## RELATED MODELS

(models in the NetLogo Models Library and elsewhere which are of related interest)

## CREDITS AND REFERENCES

(a reference to the model's URL on the web if it has one, as well as any other necessary credits, citations, and links)
@#$#@#$#@
default
true
0
Polygon -7500403 true true 150 5 40 250 150 205 260 250

airplane
true
0
Polygon -7500403 true true 150 0 135 15 120 60 120 105 15 165 15 195 120 180 135 240 105 270 120 285 150 270 180 285 210 270 165 240 180 180 285 195 285 165 180 105 180 60 165 15

arrow
true
0
Polygon -7500403 true true 150 0 0 150 105 150 105 293 195 293 195 150 300 150

box
false
0
Polygon -7500403 true true 150 285 285 225 285 75 150 135
Polygon -7500403 true true 150 135 15 75 150 15 285 75
Polygon -7500403 true true 15 75 15 225 150 285 150 135
Line -16777216 false 150 285 150 135
Line -16777216 false 150 135 15 75
Line -16777216 false 150 135 285 75

bug
true
0
Circle -7500403 true true 96 182 108
Circle -7500403 true true 110 127 80
Circle -7500403 true true 110 75 80
Line -7500403 true 150 100 80 30
Line -7500403 true 150 100 220 30

butterfly
true
0
Polygon -7500403 true true 150 165 209 199 225 225 225 255 195 270 165 255 150 240
Polygon -7500403 true true 150 165 89 198 75 225 75 255 105 270 135 255 150 240
Polygon -7500403 true true 139 148 100 105 55 90 25 90 10 105 10 135 25 180 40 195 85 194 139 163
Polygon -7500403 true true 162 150 200 105 245 90 275 90 290 105 290 135 275 180 260 195 215 195 162 165
Polygon -16777216 true false 150 255 135 225 120 150 135 120 150 105 165 120 180 150 165 225
Circle -16777216 true false 135 90 30
Line -16777216 false 150 105 195 60
Line -16777216 false 150 105 105 60

car
false
0
Polygon -7500403 true true 300 180 279 164 261 144 240 135 226 132 213 106 203 84 185 63 159 50 135 50 75 60 0 150 0 165 0 225 300 225 300 180
Circle -16777216 true false 180 180 90
Circle -16777216 true false 30 180 90
Polygon -16777216 true false 162 80 132 78 134 135 209 135 194 105 189 96 180 89
Circle -7500403 true true 47 195 58
Circle -7500403 true true 195 195 58

circle
false
0
Circle -7500403 true true 0 0 300

circle 2
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240

cow
false
0
Polygon -7500403 true true 200 193 197 249 179 249 177 196 166 187 140 189 93 191 78 179 72 211 49 209 48 181 37 149 25 120 25 89 45 72 103 84 179 75 198 76 252 64 272 81 293 103 285 121 255 121 242 118 224 167
Polygon -7500403 true true 73 210 86 251 62 249 48 208
Polygon -7500403 true true 25 114 16 195 9 204 23 213 25 200 39 123

cylinder
false
0
Circle -7500403 true true 0 0 300

dot
false
0
Circle -7500403 true true 90 90 120

face happy
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 255 90 239 62 213 47 191 67 179 90 203 109 218 150 225 192 218 210 203 227 181 251 194 236 217 212 240

face neutral
false
0
Circle -7500403 true true 8 7 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Rectangle -16777216 true false 60 195 240 225

face sad
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 168 90 184 62 210 47 232 67 244 90 220 109 205 150 198 192 205 210 220 227 242 251 229 236 206 212 183

fish
false
0
Polygon -1 true false 44 131 21 87 15 86 0 120 15 150 0 180 13 214 20 212 45 166
Polygon -1 true false 135 195 119 235 95 218 76 210 46 204 60 165
Polygon -1 true false 75 45 83 77 71 103 86 114 166 78 135 60
Polygon -7500403 true true 30 136 151 77 226 81 280 119 292 146 292 160 287 170 270 195 195 210 151 212 30 166
Circle -16777216 true false 215 106 30

flag
false
0
Rectangle -7500403 true true 60 15 75 300
Polygon -7500403 true true 90 150 270 90 90 30
Line -7500403 true 75 135 90 135
Line -7500403 true 75 45 90 45

flower
false
0
Polygon -10899396 true false 135 120 165 165 180 210 180 240 150 300 165 300 195 240 195 195 165 135
Circle -7500403 true true 85 132 38
Circle -7500403 true true 130 147 38
Circle -7500403 true true 192 85 38
Circle -7500403 true true 85 40 38
Circle -7500403 true true 177 40 38
Circle -7500403 true true 177 132 38
Circle -7500403 true true 70 85 38
Circle -7500403 true true 130 25 38
Circle -7500403 true true 96 51 108
Circle -16777216 true false 113 68 74
Polygon -10899396 true false 189 233 219 188 249 173 279 188 234 218
Polygon -10899396 true false 180 255 150 210 105 210 75 240 135 240

house
false
0
Rectangle -7500403 true true 45 120 255 285
Rectangle -16777216 true false 120 210 180 285
Polygon -7500403 true true 15 120 150 15 285 120
Line -16777216 false 30 120 270 120

leaf
false
0
Polygon -7500403 true true 150 210 135 195 120 210 60 210 30 195 60 180 60 165 15 135 30 120 15 105 40 104 45 90 60 90 90 105 105 120 120 120 105 60 120 60 135 30 150 15 165 30 180 60 195 60 180 120 195 120 210 105 240 90 255 90 263 104 285 105 270 120 285 135 240 165 240 180 270 195 240 210 180 210 165 195
Polygon -7500403 true true 135 195 135 240 120 255 105 255 105 285 135 285 165 240 165 195

line
true
0
Line -7500403 true 150 0 150 300

line half
true
0
Line -7500403 true 150 0 150 150

pentagon
false
0
Polygon -7500403 true true 150 15 15 120 60 285 240 285 285 120

person
false
0
Circle -7500403 true true 110 5 80
Polygon -7500403 true true 105 90 120 195 90 285 105 300 135 300 150 225 165 300 195 300 210 285 180 195 195 90
Rectangle -7500403 true true 127 79 172 94
Polygon -7500403 true true 195 90 240 150 225 180 165 105
Polygon -7500403 true true 105 90 60 150 75 180 135 105

plant
false
0
Rectangle -7500403 true true 135 90 165 300
Polygon -7500403 true true 135 255 90 210 45 195 75 255 135 285
Polygon -7500403 true true 165 255 210 210 255 195 225 255 165 285
Polygon -7500403 true true 135 180 90 135 45 120 75 180 135 210
Polygon -7500403 true true 165 180 165 210 225 180 255 120 210 135
Polygon -7500403 true true 135 105 90 60 45 45 75 105 135 135
Polygon -7500403 true true 165 105 165 135 225 105 255 45 210 60
Polygon -7500403 true true 135 90 120 45 150 15 180 45 165 90

sheep
false
15
Circle -1 true true 203 65 88
Circle -1 true true 70 65 162
Circle -1 true true 150 105 120
Polygon -7500403 true false 218 120 240 165 255 165 278 120
Circle -7500403 true false 214 72 67
Rectangle -1 true true 164 223 179 298
Polygon -1 true true 45 285 30 285 30 240 15 195 45 210
Circle -1 true true 3 83 150
Rectangle -1 true true 65 221 80 296
Polygon -1 true true 195 285 210 285 210 240 240 210 195 210
Polygon -7500403 true false 276 85 285 105 302 99 294 83
Polygon -7500403 true false 219 85 210 105 193 99 201 83

square
false
0
Rectangle -7500403 true true 30 30 270 270

square 2
false
0
Rectangle -7500403 true true 30 30 270 270
Rectangle -16777216 true false 60 60 240 240

star
false
0
Polygon -7500403 true true 151 1 185 108 298 108 207 175 242 282 151 216 59 282 94 175 3 108 116 108

target
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240
Circle -7500403 true true 60 60 180
Circle -16777216 true false 90 90 120
Circle -7500403 true true 120 120 60

tree
false
0
Circle -7500403 true true 118 3 94
Rectangle -6459832 true false 120 195 180 300
Circle -7500403 true true 65 21 108
Circle -7500403 true true 116 41 127
Circle -7500403 true true 45 90 120
Circle -7500403 true true 104 74 152

triangle
false
0
Polygon -7500403 true true 150 30 15 255 285 255

triangle 2
false
0
Polygon -7500403 true true 150 30 15 255 285 255
Polygon -16777216 true false 151 99 225 223 75 224

truck
false
0
Rectangle -7500403 true true 4 45 195 187
Polygon -7500403 true true 296 193 296 150 259 134 244 104 208 104 207 194
Rectangle -1 true false 195 60 195 105
Polygon -16777216 true false 238 112 252 141 219 141 218 112
Circle -16777216 true false 234 174 42
Rectangle -7500403 true true 181 185 214 194
Circle -16777216 true false 144 174 42
Circle -16777216 true false 24 174 42
Circle -7500403 false true 24 174 42
Circle -7500403 false true 144 174 42
Circle -7500403 false true 234 174 42

turtle
true
0
Polygon -10899396 true false 215 204 240 233 246 254 228 266 215 252 193 210
Polygon -10899396 true false 195 90 225 75 245 75 260 89 269 108 261 124 240 105 225 105 210 105
Polygon -10899396 true false 105 90 75 75 55 75 40 89 31 108 39 124 60 105 75 105 90 105
Polygon -10899396 true false 132 85 134 64 107 51 108 17 150 2 192 18 192 52 169 65 172 87
Polygon -10899396 true false 85 204 60 233 54 254 72 266 85 252 107 210
Polygon -7500403 true true 119 75 179 75 209 101 224 135 220 225 175 261 128 261 81 224 74 135 88 99

wheel
false
0
Circle -7500403 true true 3 3 294
Circle -16777216 true false 30 30 240
Line -7500403 true 150 285 150 15
Line -7500403 true 15 150 285 150
Circle -7500403 true true 120 120 60
Line -7500403 true 216 40 79 269
Line -7500403 true 40 84 269 221
Line -7500403 true 40 216 269 79
Line -7500403 true 84 40 221 269

wolf
false
0
Polygon -16777216 true false 253 133 245 131 245 133
Polygon -7500403 true true 2 194 13 197 30 191 38 193 38 205 20 226 20 257 27 265 38 266 40 260 31 253 31 230 60 206 68 198 75 209 66 228 65 243 82 261 84 268 100 267 103 261 77 239 79 231 100 207 98 196 119 201 143 202 160 195 166 210 172 213 173 238 167 251 160 248 154 265 169 264 178 247 186 240 198 260 200 271 217 271 219 262 207 258 195 230 192 198 210 184 227 164 242 144 259 145 284 151 277 141 293 140 299 134 297 127 273 119 270 105
Polygon -7500403 true true -1 195 14 180 36 166 40 153 53 140 82 131 134 133 159 126 188 115 227 108 236 102 238 98 268 86 269 92 281 87 269 103 269 113

x
false
0
Polygon -7500403 true true 270 75 225 30 30 225 75 270
Polygon -7500403 true true 30 75 75 30 270 225 225 270
@#$#@#$#@
NetLogo 6.1.0
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
<experiments>
  <experiment name="prediction-day10-varying-mono-init" repetitions="3" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="240"/>
    <metric>getSeed</metric>
    <metric>getViability</metric>
    <metric>getRemainingCellRatio</metric>
    <enumeratedValueSet variable="gui-prop-mono-init">
      <value value="0.15"/>
      <value value="0.16"/>
      <value value="0.22"/>
      <value value="0.26"/>
      <value value="0.52"/>
      <value value="0.53"/>
      <value value="0.55"/>
      <value value="0.6"/>
      <value value="0.65"/>
      <value value="0.75"/>
      <value value="0.8"/>
      <value value="0.86"/>
      <value value="0.9"/>
      <value value="0.9"/>
      <value value="1.04"/>
      <value value="1.1"/>
      <value value="1.16"/>
      <value value="1.2"/>
      <value value="1.21"/>
      <value value="1.6"/>
      <value value="1.84"/>
      <value value="2.15"/>
      <value value="2.5"/>
      <value value="2.85"/>
      <value value="3.12"/>
      <value value="3.33"/>
      <value value="4.26"/>
      <value value="5.72"/>
      <value value="6.57"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="prediction-patient1-1.1" repetitions="8" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="312"/>
    <metric>getSeed</metric>
    <metric>getViability</metric>
    <metric>getRemainingCellRatio</metric>
    <enumeratedValueSet variable="gui-prop-mono-init">
      <value value="1.1"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="prediction-patient2-2.5" repetitions="8" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="312"/>
    <metric>getSeed</metric>
    <metric>getViability</metric>
    <metric>getRemainingCellRatio</metric>
    <enumeratedValueSet variable="gui-prop-mono-init">
      <value value="2.5"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="prediction-patient3-1.25" repetitions="8" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="312"/>
    <metric>getSeed</metric>
    <metric>getViability</metric>
    <metric>getRemainingCellRatio</metric>
    <enumeratedValueSet variable="gui-prop-mono-init">
      <value value="2.5"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="Explo-stochasticity-100-simulations" repetitions="12" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="336"/>
    <metric>getSeed</metric>
    <metric>getViability</metric>
    <metric>getRemainingCellRatio</metric>
    <enumeratedValueSet variable="gui-apo-mov">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="gui-need-sig-mov">
      <value value="9"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="gui-layers">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="gui-alpha">
      <value value="252"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="gui-mono-phago-eff">
      <value value="51"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="gui-NLC-phago-eff">
      <value value="46"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="gui-M-phago-eff">
      <value value="78"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="gui-M-kill-eff">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="gui-cll-sens-dist">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="gui-mono-sens-dist">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="gui-nlc-sens-dist">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="gui-macro-sens-dist">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="gui-nlc-threshold">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="gui-sig-init">
      <value value="70"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="gui-sig-init-std">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="gui-diff-mean">
      <value value="36"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="gui-diff-std">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="gui-life-init-gamma">
      <value value="2094"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="gui-alpha-distrib">
      <value value="0.3318906134658388"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="PerturbM2PhagoEff2" repetitions="3" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="312"/>
    <metric>getSeed</metric>
    <metric>getViability</metric>
    <metric>getRemainingCellRatio</metric>
    <steppedValueSet variable="gui-M2-phago-eff" first="0" step="2" last="10"/>
  </experiment>
  <experiment name="PerturbMonoPhagoEff" repetitions="3" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="312"/>
    <metric>getSeed</metric>
    <metric>getViability</metric>
    <metric>getRemainingCellRatio</metric>
    <steppedValueSet variable="gui-mono-phago-eff" first="0" step="2" last="10"/>
  </experiment>
  <experiment name="PerturbNLCPhagoEff" repetitions="3" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="312"/>
    <metric>getSeed</metric>
    <metric>getViability</metric>
    <metric>getRemainingCellRatio</metric>
    <steppedValueSet variable="gui-nlc-phago-eff" first="0" step="2" last="10"/>
  </experiment>
  <experiment name="PerturbM2KillEff" repetitions="3" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="312"/>
    <metric>getSeed</metric>
    <metric>getViability</metric>
    <metric>getRemainingCellRatio</metric>
    <steppedValueSet variable="gui-M2-kill-eff" first="0" step="1" last="5"/>
  </experiment>
  <experiment name="PerturbNLCProtection" repetitions="3" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="312"/>
    <metric>getSeed</metric>
    <metric>getViability</metric>
    <metric>getRemainingCellRatio</metric>
    <steppedValueSet variable="gui-alpha" first="0" step="10" last="50"/>
  </experiment>
  <experiment name="PerturbApoMov" repetitions="3" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="312"/>
    <metric>getSeed</metric>
    <metric>getViability</metric>
    <metric>getRemainingCellRatio</metric>
    <steppedValueSet variable="gui-apo-mov" first="0" step="2" last="10"/>
  </experiment>
  <experiment name="PerturbNeedSigMov" repetitions="3" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="312"/>
    <metric>getSeed</metric>
    <metric>getViability</metric>
    <metric>getRemainingCellRatio</metric>
    <steppedValueSet variable="gui-need-sig-mov" first="0" step="2" last="10"/>
  </experiment>
  <experiment name="PerturbDiffStd" repetitions="3" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="312"/>
    <metric>getSeed</metric>
    <metric>getViability</metric>
    <metric>getRemainingCellRatio</metric>
    <steppedValueSet variable="gui-diff-std" first="0" step="10" last="50"/>
  </experiment>
  <experiment name="PerturbDiffMean" repetitions="3" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="312"/>
    <metric>getSeed</metric>
    <metric>getViability</metric>
    <metric>getRemainingCellRatio</metric>
    <steppedValueSet variable="gui-diff-mean" first="0" step="24" last="120"/>
  </experiment>
  <experiment name="PerturbLayers" repetitions="3" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="312"/>
    <metric>getSeed</metric>
    <metric>getViability</metric>
    <metric>getRemainingCellRatio</metric>
    <steppedValueSet variable="gui-layers" first="1" step="1" last="2"/>
  </experiment>
  <experiment name="PerturbCLLSensDist" repetitions="3" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="312"/>
    <metric>getSeed</metric>
    <metric>getViability</metric>
    <metric>getRemainingCellRatio</metric>
    <steppedValueSet variable="gui-cll-sens-dist" first="1" step="1" last="3"/>
  </experiment>
  <experiment name="PerturbMonoSensDist" repetitions="3" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="312"/>
    <metric>getSeed</metric>
    <metric>getViability</metric>
    <metric>getRemainingCellRatio</metric>
    <steppedValueSet variable="gui-mono-sens-dist" first="1" step="1" last="3"/>
  </experiment>
  <experiment name="PerturbNLCSensDist" repetitions="3" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="312"/>
    <metric>getSeed</metric>
    <metric>getViability</metric>
    <metric>getRemainingCellRatio</metric>
    <steppedValueSet variable="gui-nlc-sens-dist" first="1" step="1" last="3"/>
  </experiment>
  <experiment name="PerturbMacroSensDist" repetitions="3" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="312"/>
    <metric>getSeed</metric>
    <metric>getViability</metric>
    <metric>getRemainingCellRatio</metric>
    <steppedValueSet variable="gui-macro-sens-dist" first="1" step="1" last="3"/>
  </experiment>
  <experiment name="PerturbNLCThreshold" repetitions="3" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="312"/>
    <metric>getSeed</metric>
    <metric>getViability</metric>
    <metric>getRemainingCellRatio</metric>
    <steppedValueSet variable="gui-nlc-threshold" first="0" step="30" last="180"/>
  </experiment>
  <experiment name="PerturbSigInit" repetitions="3" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="312"/>
    <metric>getSeed</metric>
    <metric>getViability</metric>
    <metric>getRemainingCellRatio</metric>
    <steppedValueSet variable="gui-sig-init" first="0" step="30" last="120"/>
  </experiment>
  <experiment name="PerturbSigInitStd" repetitions="3" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="312"/>
    <metric>getSeed</metric>
    <metric>getViability</metric>
    <metric>getRemainingCellRatio</metric>
    <steppedValueSet variable="gui-sig-init-std" first="0" step="10" last="50"/>
  </experiment>
  <experiment name="prediction-day8-varying-mono-init" repetitions="3" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="192"/>
    <metric>getSeed</metric>
    <metric>getViability</metric>
    <metric>getRemainingCellRatio</metric>
    <enumeratedValueSet variable="gui-prop-mono-init">
      <value value="0.8"/>
      <value value="0.26"/>
      <value value="0.9"/>
      <value value="6.57"/>
      <value value="1.6"/>
      <value value="0.15"/>
      <value value="0.16"/>
      <value value="1.21"/>
      <value value="1.16"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="prediction-day12-varying-mono-init" repetitions="3" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="288"/>
    <metric>getSeed</metric>
    <metric>getViability</metric>
    <metric>getRemainingCellRatio</metric>
    <enumeratedValueSet variable="gui-prop-mono-init">
      <value value="0.8"/>
      <value value="0.26"/>
      <value value="0.9"/>
      <value value="6.57"/>
      <value value="1.6"/>
      <value value="0.15"/>
      <value value="0.16"/>
      <value value="1.21"/>
      <value value="1.16"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="prediction-day9-varying-mono-init" repetitions="3" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="216"/>
    <metric>getSeed</metric>
    <metric>getViability</metric>
    <metric>getRemainingCellRatio</metric>
    <enumeratedValueSet variable="gui-prop-mono-init">
      <value value="0"/>
      <value value="0.16"/>
      <value value="0.31"/>
      <value value="0.63"/>
      <value value="1.25"/>
      <value value="2.5"/>
      <value value="5"/>
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="gui-apo-mov">
      <value value="8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="gui-need-sig-mov">
      <value value="9"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="gui-layers">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="gui-alpha">
      <value value="244"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="gui-mono-phago-eff">
      <value value="25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="gui-NLC-phago-eff">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="gui-M-phago-eff">
      <value value="92"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="gui-M-kill-eff">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="gui-cll-sens-dist">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="gui-mono-sens-dist">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="gui-nlc-sens-dist">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="gui-macro-sens-dist">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="gui-nlc-threshold">
      <value value="99"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="gui-sig-init">
      <value value="63"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="gui-sig-init-std">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="gui-diff-mean">
      <value value="49"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="gui-diff-std">
      <value value="6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="gui-life-init-gamma">
      <value value="2120"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="gui-alpha-distrib">
      <value value="0.3188403442070939"/>
    </enumeratedValueSet>
  </experiment>
</experiments>
@#$#@#$#@
@#$#@#$#@
default
0.0
-0.2 0 0.0 1.0
0.0 1 1.0 0.0
0.2 0 0.0 1.0
link direction
true
0
Line -7500403 true 150 150 90 180
Line -7500403 true 150 150 210 180
@#$#@#$#@
0
@#$#@#$#@
