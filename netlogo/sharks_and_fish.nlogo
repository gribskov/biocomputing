breed [ fishes fish ]
breed [ sharks shark ]

fishes-own [ love hate backstab panic-amt energy age ]

; globals defined in interface
;    fish_n       number of fishes to create
;    view_limit   how far the fishes and shark can see

globals [ fishcount maxpanic mutation_rate energy_increment fishspeed sharkspeed]
;     fishcount         total number of fish created in simulation
;     maxpanic          maximum distance a fish moves when panicing
;     mutation_rate     probability that a new fish is mutated
;     energy_increment  amount of energy fish gain at each tick
;     fishspeed         speed of the fishes
;     sharkspeed        speed of the shark

;--------------------------------------------------------------------------------
; go button
;--------------------------------------------------------------------------------
to set-parameters
  set maxpanic  100
  set mutation_rate 0.1
  set energy_increment 0.05
  set fishspeed 0.8
  set sharkspeed 1.0
end

;--------------------------------------------------------------------------------
; setup button
;--------------------------------------------------------------------------------
to setup
  clear-all
  set-parameters

  create-fishes fish_n [
    make-fish
  ]
  set fishcount fish_n

  create-sharks 1 [
    setxy random-xcor random-ycor
    set color red
    set shape "arrow"
    set size 6
    ;pen-down
  ]
  reset-ticks
  tick

end

;--------------------------------------------------------------------------------
; main event loop
;--------------------------------------------------------------------------------
to go

  move-sharks
  move-fishes
  panic
  eat-fish
  tick

end

;--------------------------------------------------------------------------------
; find the closest target (fish or shark)  and if it is less than view_limit
; distance away, return the dx, dy, and the id of the target
; if it's outside the distance return 0, 0, target_id
;--------------------------------------------------------------------------------
to-report min-distance [ target scale ]
  let closest-target min-one-of target [ distance self]
  let d distance closest-target
;  show d
  ifelse d > view_limit
    [ report ( list 0 0 any? target ) ]
    [ let ddx ([xcor] of closest-target) - xcor
      let ddy ([ycor] of closest-target) - ycor
      ;show list (xcor) (ycor)
;      show (list (view_limit) (ddx) (ddy))
      report (list (ddx * scale) (ddy * scale) closest-target ) ]
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; shark procedures
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;--------------------------------------------------------------------------------
; the shark moves toward the closest fish within its view_limit, or turns in a
; random direction in the range -5 to 5 degrees
;--------------------------------------------------------------------------------
to move-sharks
  ask sharks [
    ifelse distance min-one-of fishes [ distance myself ] < view_limit
        [face min-one-of fishes [ distance myself]]
        [rt (random 10) - 5]
    fd sharkspeed
  ]
end

;--------------------------------------------------------------------------------
; when the fish is within 2.5, it is killed and a new one created
;--------------------------------------------------------------------------------
to eat-fish
  ask sharks [
    let prey ( min-one-of fishes [distance myself] )
    if (distance  prey ) < 2.5 [
      ask prey [ die ]
      breed-fish
    ]
  ]
end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; fish procedures
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;--------------------------------------------------------------------------------
; construct a single fish
;--------------------------------------------------------------------------------
to make-fish
    set shape "circle"
    setxy random-pxcor random-pycor
    facexy random-xcor random-ycor
    set color blue
    set size 4
    set love ( random-float 1.0 ) - 0.5
    set hate ( random-float 1.0 )
    set backstab ( random-float 180 )
    set panic-amt 0.25
    set energy 10

    set color color-fish
end

;--------------------------------------------------------------------------------
; move the fish depending on location of the the closest fish
; the closest shark (usually just one), and each fishes love and hate
; genes. if neither fish nor shark is in sight, continue forward with
; 0 - 10 degree random turn]
;--------------------------------------------------------------------------------
to move-fishes
  ask fishes [
    let closest-fish min-distance fishes love
    let closest-shark min-distance sharks hate
    let xnew item 0 closest-fish - item 0 closest-shark
    let ynew item 1 closest-fish - item 1 closest-shark
    facexy (xcor + xnew) (ycor + ynew)

    lt ((random 20) - 10 )
    fd fishspeed

    set energy ( energy + energy_increment )
    if energy > maxpanic [ set energy maxpanic ]
    set age ( age + 1 )
  ]
end

;--------------------------------------------------------------------------------
; when the shark gets close the fish can use a burst of energy to try escaping
; panic-amt is the fraction of energy to use, but the amount moved cannot be
; more than maxpanic
; the fish will move toward the closest fish, offset by backstab degrees to the
; right
;--------------------------------------------------------------------------------
to panic
  ask fishes [
    if distance min-one-of sharks [ distance myself ] < 20 [
      if energy > 0 [
        let panicdist ( energy * panic-amt )
        if panicdist > maxpanic [ set panicdist maxpanic ]
        face min-one-of other fishes [ distance myself ]
        ;rt ( random backstab - ( backstab / 2 ) )
        rt backstab
        fd panicdist
        set energy (energy - panicdist )
      ]
    ]
  ]
end

;--------------------------------------------------------------------------------
; randomly choose an existing fish and duplicate
; mutation_rate is the probability of introducting mutations
; if mutations are made, love, hate, backstab, and panic-amt are altered
; by a random amount betweeen -10 and +10 %
;--------------------------------------------------------------------------------
to breed-fish
  ask one-of fishes [
    set fishcount ( fishcount + 1 )
    hatch 1 [
      set age 0
      setxy random-xcor random-ycor
      if mutation_rate > random-float 1.0 [
        set love ( love + random-float 0.2 - 0.1 )
        if love > 1.0 [ set love 1.0 ]
        if love < -1.0 [ set love -1.0 ]

        set hate ( hate + random-float 0.2 - 0.1 )
        if hate > 1.0 [ set hate 1.0 ]
        if hate < 0.0 [ set hate 0.0 ]

        set backstab ( backstab + random-float 20 - 10 )
        if backstab > 360 [ set backstab 360 ]
        if backstab < 0 [ set backstab 0 ]

        set panic-amt ( panic-amt + random-float 0.2 - 0.1 )
        if panic-amt > 1.0 [ set panic-amt 1.0 ]
        if panic-amt < 0.0 [ set panic-amt 0.0 ]

      ]
      set color color-fish
      fd fishspeed
    ]
  ]
end

;--------------------------------------------------------------------------------
; determine color based on fish love, hate, and panic-amt.
; reports an RGB color
;--------------------------------------------------------------------------------
to-report color-fish
  let r 255 * (love + 1) / 2
  let b 255 * hate
  let g 255 * panic-amt
  report rgb r g b
end
