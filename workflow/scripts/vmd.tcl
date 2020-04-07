set infile [lindex $argv 0]
set gifdir [lindex $argv 1]

mol new ${infile} type {xyz} first 0 last -1 step 1 waitfor 1
mol modcolor 0 0 Element
mol modstyle 0 0 VDW 1.000000 12.000000
mol modmaterial 0 0 Opaque

file mkdir ${gifdir}

set numframes [molinfo top get numframes]
set frame 0

# Set axis spin per timestep
set x_inc 1
set y_inc 1
set z_inc 1

for {set i 0} {$i < $numframes} {incr i} {
  animate goto $i
  set filename ${gifdir}/snap.[format "%04d" $frame].tga
  render TachyonInternal ${filename}
  incr frame
  rotate x by ${x_inc}
  rotate y by ${y_inc}
  rotate z by ${z_inc}
}

set angle_inc_end 4
for {set i 0} {$i < 360} {incr i ${angle_inc_end}} {
  set filename ${gifdir}/snap.[format "%04d" $frame].tga
  render TachyonInternal ${filename}
  incr frame
  rotate y by ${angle_inc_end}
}
