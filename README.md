# space-time-globe
An interactive space-time globe written in Python

The space-time globe is an idea that I initially saw on minutephysics: https://www.youtube.com/watch?v=Rh0pYtQG5wI
This globe is such an intuitive and downright revolutionary way to present special relativity. The Lorentz transformation
comes to life and is far easier to understand and play with than using pure maths.

This initial version is not refined in terms of user interaction, but does the trick of showing the basics. It draws a 
space-time diagram with the time cone and a set of proper time curves for s=1,2,3,4,... A number of events defined in 
the program are shown with different colors. Each points simultaneity lines are shown with dotted lines in the same color
as the event. The legend shows the events and their length in proper time.

There is only one python file: SpaceTimeGlobe.py, which will run and present the graphics (using python3).
Events are entered directly in the program near the end where the 'events = ...' statement is located.
The mouse is used to drag the whole scenario around and change the reference frame.
The keys '1'-'3' may be used to change directly to the perspective of events 1-3.
The "'" key may be used to toggle on/off the grid (mesh).
The ";" key toogles on/off the simultaneity lines.
The "." key toggles on/off the proper time lines.

Thats it for now. The default scenario shows a variant of the twin paradox.
