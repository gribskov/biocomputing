patches-own [ nbors ]

to setup
  clear-all
  set evolve false
  reset-ticks
end

to go
  check-mouse
  if evolve [
    count-neighbors
    color-patches
    update-cells
    tick
  ]
end

to check-mouse
  if mouse-down? [
    ask patch round mouse-xcor round mouse-ycor [ set pcolor red ]
    display
  ]
end

to count-neighbors
  ask patches [
    set nbors count neighbors with [ pcolor = red ]
  ]
end

to color-patches
  ask patches [
    ifelse pcolor = red
      [ ]
      [ set pcolor scale-color green nbors 0 8 ]
  ]
end

to update-cells
  ask patches [
    ifelse nbors > 1 and nbors < 4
      [ ifelse pcolor = red
        [ ]
        [ if nbors > 2 [ set pcolor red ] ]
      ]
      [set pcolor scale-color green nbors 0 8 ]
    ]

end
