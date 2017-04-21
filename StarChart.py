# StarChart -- Sky chart program for visual astronomy
#
# FIXME: error in time or zone is automatically corrected without notification.
# FIXME: the notify alerts have no callback.  (Fixing this could be helpful in fixing the above.)
#
# Copyright (c) 2008 by David A. Wallace
# Released by the author under GPL license
#
# Version 0.54 (beta 5) of 2008.09.24.2210 UT
#
# See http://wiki.laptop.org/go/StarChart for operating instructions.
#
# The catalog data and algorithms herein were gleaned from several sources:
# Practical Astronomy with Your Calculator
#   by Peter Duffett-Smith (3rd ed.)
#
# Keplerian Elements for Approximate Positions of the Major Planets
#   http://ssd.jpl.nasa.gov/txt/aprx_pos_planets.pdf
#
# Bright Star Catalog
#   http://adc.gsfc.nasa.gov/adc-cgi/cat.pl?/catalogs/5/5050/
#
# Names of Stars
#   extracted from data at http://www.alcyone.de/SIT/bsc/bsc.html
#
# Messier Catalog
#   extracted from data at http://astro.nineplanets.org
#
# Constellation figures -- derived from charts at
#   http://www.astro.wisc.edu/~dolan/constellations/
# and the coordinates of the stars that the line-segments interconnect.
#
# -------------------------------------------------------------------------------
#
# INDEX
#
# (Line numbers given below are approximate, within 25 lines, usually.)
#
#   Line No.    Content
#   --------    --------------------------------------
#      100      Planetary Catalog data
#      150      Start of the code -- astronomical algorithms
#      550      Definitions for the GUI controls (toolbars)
#      625      Definition for the ChartDisplay (main GUI) class
#     1575      Definition for the StarChartActivity class
#


# ================== NON-STANDARD BUT REQUIRED LOCALIZATION =====================


# The settings for the observer's coordinates
try:
  import observatory
  longitude = observatory.data[0]
  latitude = observatory.data[1]
except:
  longitude = 0.0
  latitude = 0.0

# "nonlocal_timezone_correction" allows your XO to keep one timezone while the
# program's "now" expresses time in another timezone.  This value is normally set
# to zero, but if you are traveling, you might find it more convenient to leave
# the XO's local timekeeping alone and simply offset the star chart's time zone
# so that "now" is correct for your current locale.
#
nonlocal_timezone_correction = 0.0

# If traveling, create "travel.py" to define alternate location and timezone.
# (If "travel.py" does not exist, we use our home observatory parameters.)

try:
  import travel
  longitude = travel.data[0]
  latitude = travel.data[1]
  nonlocal_timezone_correction = travel.data[2]
except:
  pass


# =================================== IMPORTS ===================================


import pygtk
pygtk.require('2.0')
import gtk
import sys
import os
import gobject
from datetime import *
from time import sleep
from time import time as systime
from math import sin, cos, tan, asin, acos, atan, pi, sqrt, atan2
from sugar.activity import activity
from sugar.activity.activity import get_bundle_path
from sugar.graphics.alert import Alert
from sugar.graphics.alert import NotifyAlert
import logging
from gettext import gettext

# Defensive method of gettext use
def _(s):
  istrsTest = {}
  for i in range (0,4):
    istrsTest[str(i)] = str(i)

  try:
    i = gettext(s)
  except:
    i = s
  return i


# ================================= CATALOG DATA ================================
#
# Keplarian data for planets (J2000.0 epoch coordinates, except moon is J1990.5)
#
# FIXME: add data for dwbar, dI and dO
#
#   planet     arg of peri  eccentricity radius (au) inclination    RAAN        mean longitude   delta mean longitude
#   pname         wbar0       E            a            I0            O0             L                dL
planets = [
  ('Mercuy',    77.45779628, 0.20563593,  0.38709927,  7.00497902,  48.33076593, 252.25032350,   149472.67411175),
  ('Venus',    131.60246718, 0.00677672,  0.72333566,  3.39467605,  76.67984255, 181.97909950,    58517.81538729),
  ('Earth',    102.93768193, 0.01671123,  1.00000261, -0.00001531,   0.0,        100.46457166,    35999.37244981),
  ('Mars',     -23.94362959, 0.09339410,  1.52371034,  1.84969142,  49.55953891,  -4.55343205,    19140.30268499),
  ('Jupiter',   14.72847983, 0.04838624,  5.20288700,  1.30439695, 100.47390909,  34.39644051,     3034.74612775),
  ('Saturn',    92.59887831, 0.05386179,  9.53667594,  2.48599187, 113.66242448,  49.95424423,     1222.49362201),
  ('Uranus',   170.95427630, 0.04725744, 19.18916464,  0.77263783,  74.01692503, 313.23810451,      428.48202785),
  ]
sun = ('Sun',  282.93768193, 0.01671123,  1.00000261,  0.0,          0.0,        280.46457166,    35999.37244981)
#       mean longitude lon of peri   lon of Node   Inclination    Eccentricity   radius (km)   Parallax   Offset of Epoch
#           L0             P0            N0             I             e            a              phi0        tau
moon = ('Moon', 318.351648,   36.340410,      318.510107,      5.145396,   0.054900,     384401,          0.9507,    2447891.5)

# obliquity of J2000 epoch is 23.43928 degrees -- we need this value (in radians) for planetary calculations
eps = pi * 23.43928 / 180.0

# The star catalog is imported from stars1.py.
# (I'm hoping to be able to import supplementary catalogs and just add them.)

import stars1
star_chart = stars1.data
import constellations
figures = constellations.data


# -------------------------------------------------------------------------------

try:
  import dso1
  dso_chart = dso1.data
except:
  dso_chart = []



# ============================= STATE INFORMATION ===============================

# settings for display options
nightvision = False
invertdisplay = False
fliphorizontally = False
drawconstellations = True
limitingmagnitude = 4.0
# settings for time
specifytime = False
now = datetime.utcnow()
zoneoffset = -300


# ============================ START OF CODE ====================================


# Because Python's trig is done in radians and I often need answers in degrees,
# I often need to convert between radians and degrees.  These two functions are
# for doing that.

def dtor(a):
  return a * pi / 180.0

def rtod(r):
  return r * 180.0 / pi


# This function converts a decimal time to hours, minutes and seconds format and
# returns a time object.

def floattotime(f):
  h = int(f)
  mm = (f - h) * 60.0
  m = int(mm)
  s = (mm - m) * 60.0
  s = int(s)
  t = time(h, m, s)
  return t


# This function converts a decimal angle to degrees, minutes and seconds format
# and returns a string "dddmmss".

def floattoangle(f):
  h = int(f)
  mm = (f - h) * 60.0
  m = int(mm)
  s = (mm - m) * 60.0
  s = int(s)
  return '%(degree)03dd%(minute)02dm%(second)02ds' % \
         { 'degree' : h, 'minute' : abs(m), 'second': abs(s)}


# Convert a degrees-minutes-seconds angle to fractional degrees.  (This function
# will accept any character as separator but you must specify degrees, minutes
# and seconds, e.g.: 012d34m56s, or degrees and minutes, e.g.: 012d34m.)

def angletofloat(s):
  try:
    d = 0.0
    i = 0
    while ((i < 4) and (i < len(s))):
      c = s[i]
      c = c[0]
      if (c.isdigit()):
        i = i + 1
      else:
        break
    t = s[0:i]
    d = float(t)
    i = i + 1
    if (i < len(s)):
      t = s[i:]
      s = t
      i = 0
      while ((i < 3) and (i < len(s))):
        c = s[i]
        c = c[0]
        if (c.isdigit()):
          i = i + 1
        else:
          break
      t = s[0:i]
      d = d + float(t)/60.0    
      i = i + 1
      if (i < len(s)):
        t = s[i:]
        s = t
        i = 0
        while ((i < 3) and (i < len(s))):
          c = s[i]
          c = c[0]
          if (c.isdigit()):
            i = i + 1
          else:
            break
        t = s[0:i]
        d = d + float(t)/3600.0    
    return d
  except:
    return -1.0


# Converts a utc date-time to julian day

def jtime(dd):
  y = dd.year
  m = dd.month
  d = dd.day + dd.hour / 24.0 + dd.minute / 1440.0 + dd.second / 86400.0
  if (m < 3):
    m = m + 12
    y = y - 1
  a = int(y / 100.0)
  b = 2 - a + int(a / 4)
  c = int(365.25 * y)
  cd = int(30.6001 * (m + 1))
  t = b + c + cd + d + 1720994.5
  return t


# Convert a utc date-time to gst (Greenwitch Siderial Time)

def gst(d):
  d1 = datetime(d.year, d.month, d.day, 0, 0, 0)
  t1 = d.time()
  j = jtime(d1)
  s = j - 2451545.0
  t = s / 36525.0
  t0 = 6.697374558 + 2400.051336 * t + 0.000025862 * t * t
  while (t0 < 0.0):
    t0 = t0 + 24.0
  while (t0 > 24.0):
    t0 = t0 - 24.0
  ut = t1.hour + t1.minute / 60.0 + t1.second / 3600.0
  ut = ut * 1.002737909
  t0 = t0 + ut
  while (t0 < 0.0):
    t0 = t0 + 24.0
  while (t0 > 24.0):
    t0 = t0 - 24.0
  return t0


# Convert gst to lst

def lst(g):
  l = longitude / 15.0
  lst = g + l
  while (lst < 0.0):
    lst = lst + 24.0
  while (lst > 24.0):
    lst = lst - 24.0
  return lst


# convert RA to hour angle

def hrangl(ra, ut):
  g = gst(ut)
  l = lst(g)
  h = l - ra
  while (h < 0.0):
    h = h + 24.0
  while (h > 24.0):
    h = h - 24.0
  return h


# Adjust J2000 equatorial coordinates for precession.

def epochpolartonow(polar, ut):
  dec = polar[1]
  ra = polar[0]
# this method is an approximation which fails when the object is near the
# celestial pole.  So we return the epoch's ra -- a better answer than what the
# algorithm would produce
  if (dec > 88.0):
     return (ra, dec)
  if (dec < -88.0):
     return (ra, dec)
# dec is in tolerance.
  n = (ut.year - 2000.0) + (ut.month - 1.0) / \
      12.0 + 15.0 / 360.0
  ra = ra * 15.0
  ra1 = ra + ((3.0742 + \
        1.33589 * sin(dtor(ra)) * \
        tan(dtor(dec))) / 3600.0)  * n
  dec1 = dec + (20.0383 * cos(dtor(ra))) / 3600.0 * n
  return (ra1 / 15.0, dec1)


# Convert equatorial coordinates to azimuth, altitude

def polartoazalt(polar, ut):
  ra = polar[0]
  dec = polar[1]
  h = hrangl(ra, ut)
  h = h * 15.0
  a = asin(sin(dtor(dec)) * sin(dtor(latitude)) + cos(dtor(dec)) *
           cos(dtor(latitude)) * cos(dtor(h)))
  az = acos((sin(dtor(dec)) - sin(dtor(latitude)) * sin(a)) /
            (cos(dtor(latitude)) * cos(a)))
  sh = sin(dtor(h))
  if (sh > 0.0):
    az = 2.0 * pi - az
  return (rtod(az), rtod(a))


# solve M = E - (57.29578 * e) * sin(E) for E
# (not currently used -- the planetary results are accurate enough without it.)

def eccentric_anomaly(M, e):
  dE = .05
  estar = 180.0 * e / pi
# Iterate dM = M - En - estar * sin(En); dE = dM / (1 - e * cos(En)); En = En + dE
# until abs(dE) <= 0.001
# note that M and estar are in degrees, as is E and dE so make the appropriate
# conversion to radians before doing trig.
  En = M
  while (abs(dE) > .001):
    dM = M - (En - estar * sin(dtor(En)))
    dE = dM / (1 - e * cos(dtor(En)))
    En = En + dE
  return En

  
# Determine the number of minutes from GMT specified by the zone offset string.

def parse_zone_offset(s):
  try:
    bneg = False
    oh = 0
    om = 0
    i = 0
    if (s[0] == '-'):
      bneg = True
      i = 1
    elif (s[0] == '+'):
      i = 1
    j = i
    while (s[j].isdigit()):
      j = j + 1
    if (j == i):
      return (0, 0)
    oh = int(s[i:j])
    i = j
    if (s[i] != ':'):
      return (0, 0)
    i = i + 1
    om = int(s[i:])
    if (bneg):
      oh = -oh
      om = -om
    return (oh, om)
  except:
    return (0, 0)


# Determine the datetime value represented by the timestamp string.

def parse_timestamp(s):
  try:
    if ((s[4] != '/') or (s[7] != '/') or (s[10] != ',') or (s[13] != ':')):
      return (2000, 1, 1, 0, 0)
    Y = int(s[0:4])
    M = int(s[5:7])
    D = int(s[8:10])
    h = int(s[11:13])
    m = int(s[14:])
    return (Y, M, D, h, m)
  except:
    return (2000, 1, 1, 0, 0)


# Get the GMT datetime value from the timestamp and offset strings

def get_time_and_UTC_offset(timestr, offsetstr):
  class TZ(tzinfo):
    def __init__(self, offset, name):
        self.__offset = timedelta(minutes = offset)
        self.__name = name

    def utcoffset(self, dt):
        return self.__offset

    def tzname(self, dt):
        return self.__name

    def dst(self, dt):
        return timedelta(0)

# Parse the timestamp string and the UTC offset string into year, month, day,
# hour and minute and hour and minute, respectively.

  (sy, sM, sd, sh, sm) = parse_timestamp(timestr)
  (oh, om) = parse_zone_offset(offsetstr)

# now convert the zone offset to a timezone object
             
  tzo = oh * 60 + om
  tz = TZ(tzo, "")

# using the parsed time and timezone object, construct a datetime object that
# represents the local time.
             
  lt = datetime(sy, sM, sd, sh, sm, 0, 0, tz)

# finally, subtract the zone object from the local time to get GMT.

  gt = lt - tz.utcoffset(0)
  return gt


# Using the current local time, the nonlocal_timezone_correction (if any) and
# the GMT time, set the text in entry3 and entry4 to the "real" local time and
# the "real" zone offset, respectively

def set_time_and_UTC_offset():
# TODO: convert nonlocal_timezone_correction to something that a datetime can do math with
  to = timedelta(0, 60 * nonlocal_timezone_correction)
  gt = datetime.utcnow()
  lt = datetime.fromtimestamp(systime())
  dt = lt - gt + to
  tt = lt
  tt = tt + to
  lts = tt.strftime("%Y/%m/%d,%H:%M")
  dth = dt.days * 24 + int(dt.seconds / 3600)
  dtm = (dt.seconds / 60) % 60
  utos = "%d:%02d" % (dth, dtm)
  return (lts, utos)

def syntax_check_time():
  try:
    s = entry3.get_text()
# this string must conform to YYYY/MM/DD,HH:MM format and be a valid date and time.
    if (len(s) != 16):
      return False
    if ((s[4] != '/') or (s[7] != '/') or (s[10] != ',') or (s[13] != ':')):
      return False
    z = s[0:4]
    if (not z.isdigit()):
      return False
    y = int(z)
    z = s[5:7]
    if (not z.isdigit()):
      return False
    x = int(z)
    if ((x < 1) or (x > 12)):
      return False
    z = s[8:10]
    if (not z.isdigit()):
      return False
    d = int(z)
    if ((d < 1) or (d > 31)):
      return False
    if ((x == 1) or (x == 3) or (x == 5) or (x == 7) or (x == 8) or (x == 10) or (x == 12)):
      pass
    elif (x == 2):
      if (d > 29):
        return False
      elif (((y % 4) != 0) and (d > 28)):
        return False
      elif (((y % 100) == 0) and (d > 28)):
        if ((y % 400) != 0):
          return False
    elif (d > 30):
      return False
    z = s[11:13]
    if (not z.isdigit()):
      return False
    h = int(z)
    if (h > 23):
      return False
    z = s[14]
    if (not z.isdigit()):
      return False
    m = int(z)
    if (m > 59):
      return False
    return True
  except:
    return False

def syntax_check_zone():
  try:
    s = entry4.get_text()
# this string must conform to [+|-]HH:MM and must be a valid time between -14:00 and 14:00.
    l = len(s)
    if (l < 4):
      return False
    if ((s[0] == '+') or (s[0] == '-')):
      p = 1
    else:
      p = 0
    if ((l - p) < 4):
      return False
    if (not s[p].isdigit()):
      return False
    if (s[p+1].isdigit()):
      n = 1
    else:
      n = 0
    if ((l - p - n) < 4):
      return False
    if (s[1 + p + n] != ':'):
      return False
    if ((not s[p + 2 + n].isdigit()) or (not s[p + 3 + n].isdigit())):
      return False
    z = s[p : p + 1 + n]
    if (int(z) > 14):
      return False
    z = s[p + 2 + n : p + n + 4]
    if (int(z) > 59):
      return False
    return True
  except:
    return False

# ----------------------------------------------------------------------------------------
#
# These controls affect the program's state-variables and must be global or we
# can't set or retrieve their values everywhere necessary:
#
# controls on menubar1 (_("what")):
button1 = gtk.ToggleButton(_("Night Vision"))
button2 = gtk.ToggleButton(_("Invert Display"))
button3 = gtk.ToggleButton(_("Flip L/R"))
button4 = gtk.ToggleButton(_("Draw Constellations"))
label6 = gtk.Label(_("Mag:"))
rb7 = gtk.RadioButton(None, _("0"))
rb8 = gtk.RadioButton(rb7, _("1"))
rb9 = gtk.RadioButton(rb7, _("2"))
rb10 = gtk.RadioButton(rb7, _("3"))
rb11 = gtk.RadioButton(rb7, _("4"))
rb12 = gtk.RadioButton(rb7, _("5"))
rb13 = gtk.RadioButton(rb7, _("6"))
# controls on menubar2 (_("where")):
label1 = gtk.Label(_("Longitude:"))
entry1 = gtk.Entry()
entry1.set_width_chars(10)
rb1 = gtk.RadioButton(None, _("E"))
rb2 = gtk.RadioButton(rb1, _("W"))
label2 = gtk.Label(_("Latitude:"))
entry2 = gtk.Entry()
entry2.set_width_chars(10)
rb3 = gtk.RadioButton(None, _("N"))
rb4 = gtk.RadioButton(rb3, _("S"))
button5 = gtk.Button(_("Apply location"))
# controls on menubar3 (_("when")):
rb5 = gtk.RadioButton(None, _("Now"))
rb6 = gtk.RadioButton(rb5, _("Specify"))
label4 = gtk.Label(_("Time:"))
entry3 = gtk.Entry()
entry3.set_width_chars(16)
label5 = gtk.Label(_("Offset:"))
entry4 = gtk.Entry()
entry4.set_width_chars(7)
button6 = gtk.Button(_("Apply time"))

lon_alert = Alert()
lat_alert = Alert()
time_alert = NotifyAlert()
zone_alert = NotifyAlert()

# ----------------------------------------------------------------------------------------


# Set control states and values from the state variables.

def initialize_controls():
  button1.set_active(nightvision)
  button2.set_active(invertdisplay)
  button3.set_active(fliphorizontally)
  button4.set_active(drawconstellations)
  rb13.set_active(limitingmagnitude >= 6.0)
  rb12.set_active((limitingmagnitude >= 5.0) and (limitingmagnitude < 6.0))
  rb11.set_active((limitingmagnitude >= 4.0) and (limitingmagnitude < 5.0))
  rb10.set_active((limitingmagnitude >= 3.0) and (limitingmagnitude < 4.0))
  rb9.set_active((limitingmagnitude >= 2.0) and (limitingmagnitude < 3.0))
  rb8.set_active((limitingmagnitude >= 1.0) and (limitingmagnitude < 2.0))
  rb7.set_active(limitingmagnitude < 1.0)
  entry2.set_text(floattoangle(abs(latitude)))
  rb4.set_active(latitude < 0.0)
  rb3.set_active(latitude >= 0.0)
  entry1.set_text(floattoangle(abs(longitude)))
  rb2.set_active(longitude < 0.0)
  rb1.set_active(longitude >= 0.0)
  rb6.set_active(specifytime)


# ============================== ChartDisplay Object ============================


class ChartDisplay(gtk.DrawingArea):
  def __init__(self):
    super(ChartDisplay, self).__init__()
    self.colors = {}
    self.canplot = False
    self.pangolayout = self.create_pango_layout("")
    if (not specifytime):
      gobject.timeout_add(60000, self.timer1_cb)

  def area_expose_cb(self, area, event):

# Determine the area we can draw upon and adjust the chart accordingly.

    rect = self.get_allocation()
    self.screensize = (rect.width, rect.height)
    self.margin = 40
    self.diameter = min(self.screensize[0], self.screensize[1]) - \
                    2 * self.margin
    self.xoffset = (self.screensize[0] - self.diameter) / 2 - self.margin
    self.yoffset = (self.screensize[1] - self.diameter) / 2 - self.margin

# Establish color selections (need only do this once).

    if (len(self.colors) == 0):
      self.gc = self.style.fg_gc[gtk.STATE_NORMAL]
      self.colormap = self.gc.get_colormap()
      self.colors[0] = self.colormap.alloc_color('white')
      self.colors[1] = self.colormap.alloc_color('black')
      self.colors[2] = self.colormap.alloc_color('red')
      self.colors[3] = self.colormap.alloc_color('gray')
      self.canplot = True
    self.plotchart()


  def timer1_cb(self):
    if (specifytime):
# do not redraw the chart if we're not advancing time.
      return True
    self.plotchart()
    return True

 
  def azalttoxy(self, polar):
    az = dtor(polar[0])
    alt = polar[1]

# radius is proportional to 90-alt

    r = self.diameter / 2.0 - alt * (self.diameter / 2.0 /90.0)

# convert r and az to x and y
# draw the chart so east is on the right by default.
# When flipped, east is on the left, like a map.

    if (fliphorizontally):
      x = int(self.diameter / 2.0 + r * sin(az))
    else:
      x = int(self.diameter / 2.0 - r * sin(az))
    y = int(self.diameter / 2.0 - r * cos(az))
    return (x,y)

    
  def plotchart(self):
    global now
    global zoneoffset

# Erase prior plot

    if (not self.canplot):
      return
    self.cleararea()
    if (invertdisplay):
      if (nightvision):
        self.gc.set_foreground(self.colors[2])
      else:
        self.gc.set_foreground(self.colors[0])
    else:
      self.gc.set_foreground(self.colors[1])
    self.window.draw_arc(self.gc,
                              True,
                              self.xoffset + self.margin - 2,
                              self.yoffset + self.margin - 2,
                              self.diameter + 4,
                              self.diameter + 4,
                              0,
                              360 * 64)

# Plot sky circle

    if (not invertdisplay):
      if (nightvision):
        self.gc.set_foreground(self.colors[2])
      else:
        self.gc.set_foreground(self.colors[0])
    else:
      self.gc.set_foreground(self.colors[1])
    self.window.draw_arc(self.gc,
                              False,
                              self.xoffset + self.margin - 2,
                              self.yoffset + self.margin - 2,
                              self.diameter + 4,
                              self.diameter + 4,
                              0,
                              360 * 64)

# label the cardinal points.

    if (nightvision):
      self.gc.set_foreground(self.colors[2])
    else:
      self.gc.set_foreground(self.colors[1])
    self.pangolayout.set_text(_("N"))
    self.window.draw_layout(self.gc,
                     self.xoffset + self.margin + self.diameter / 2 - 10,
                     self.margin - 30, self.pangolayout)
    self.pangolayout.set_text(_("S"))
    self.window.draw_layout(self.gc,
                     self.xoffset + self.margin + self.diameter / 2 - 10,
                     2 * self.margin + self.diameter - 30, self.pangolayout)
    if (not fliphorizontally):
      self.pangolayout.set_text(_("E"))
    else:
      self.pangolayout.set_text(_("W"))
    self.window.draw_layout(self.gc,
                     self.xoffset + self.margin - 30,
                     self.margin + self.diameter / 2 - 10, self.pangolayout)
    if (not fliphorizontally):
      self.pangolayout.set_text(_("W"))
    else:
      self.pangolayout.set_text(_("E"))
    self.window.draw_layout(self.gc,
                     self.xoffset + self.margin + self.diameter + 10,
                     self.margin + self.diameter / 2 - 10, self.pangolayout)

    if (not invertdisplay):
      if (nightvision):
        self.gc.set_foreground(self.colors[2])
      else:
        self.gc.set_foreground(self.colors[0])
    else:
      self.gc.set_foreground(self.colors[1])

# Set the time of plotting (now).

    if (not specifytime):
      now = datetime.utcnow()
      (tstr, ostr) = set_time_and_UTC_offset()
      entry3.set_text(tstr)
      entry4.set_text(ostr)
    else:
      now = get_time_and_UTC_offset(entry3.get_text(), entry4.get_text())
      (hh, mm) = parse_zone_offset(entry4.get_text())
      zoneoffset = 60 * hh
      if (hh < 0):
        zoneoffset = zoneoffset - mm
      else:
        zoneoffset = zoneoffset + mm

# Plot the stars.

    for name, (ra, dec, mag) in star_chart.iteritems():
      if (mag <= limitingmagnitude):

# convert the ra and dec from the J2000 epoch to the plot time

        polar = epochpolartonow((ra, dec), now)

# convert the equatorial coordinates to altitude and azimuth

        azalt = polartoazalt(polar, now)

# plot any object that is more than 1 degree above the horizon

        if (azalt[1] >= 0.0175):
          starsize = 2 + int(7.0 - mag)
          (px, py) = self.azalttoxy(azalt)
          px = px + self.margin - 2 + self.xoffset - starsize / 2
          py = py + self.margin - 2 + self.yoffset - starsize / 2
          self.window.draw_arc(self.gc, True,
                 px,
                 py,
                 starsize,
                 starsize,
                 0,
                 360*64)

# Plot the deep sky objects.

    for i in range(len(dso_chart)):
      (nM, strCon, ra, dec, mag, majA, minA, posA, strT, strN) = dso_chart[i]

# convert the ra and dec from the J2000 epoch to the plot time

      polar = epochpolartonow((ra, dec), now)

# convert the equatorial coordinates to altitude and azimuth

      azalt = polartoazalt(polar, now)
      if (azalt[1] >= 0.0175):
        (px, py) = self.azalttoxy(azalt)
        px = px + self.margin - 2 + self.xoffset
        py = py + self.margin - 2 + self.yoffset
        self.plot_DSO(strT, majA, minA, mag, px, py)

# Plot the constellation figures.  This is essentially the same process as for
# plotting a star but we have to figure out the alt/az coordinates for both ends
# of the line segment.

    if (drawconstellations):
      if (not invertdisplay):
        if (nightvision):
          self.gc.set_foreground(self.colors[2])
        else:
          self.gc.set_foreground(self.colors[0])
      else:
        self.gc.set_foreground(self.colors[1])
      for i in range(len(figures)):
        (ra1, dec1, ra2, dec2) = figures[i]
        polar1 = epochpolartonow((ra1, dec1), now)
        azalt1 = polartoazalt(polar1, now)
        if (azalt1[1] >=  0.0175):
          (px1, py1) = self.azalttoxy(azalt1)
          px1 = px1 + self.margin - 2 + self.xoffset
          py1 = py1 + self.margin - 2 + self.yoffset
          polar2 = epochpolartonow((ra2, dec2), now)
          azalt2 = polartoazalt(polar2, now)
          if (azalt2[1] >=  0.0175):
            (px2, py2) = self.azalttoxy(azalt2)
            px2 = px2 + self.margin - 2 + self.xoffset
            py2 = py2 + self.margin - 2 + self.yoffset
            self.window.draw_line(self.gc, px1, py1, px2, py2)

# Plot the planets, the moon and the sun.

# We need to adjust the Keplerian data by the number of centuries since the
# epoch.  Mostly, this causes the planets to revolve around the sun at the rate
# of their orbital period (N in the following computations), but the orbits
# themselves will also change over time.
# FIXME: for now, we only worry about orbital motion, not additionally orbital
#        mutation.

    T = (jtime(now) - 2451545.0) / 36525.0

# Calculate the earth's heliocentric radius and longitude coordinates so we
# can figure the geocentric ecliptic coordinates of the other planets.

    (name, wbar, e, a, I, O, L0, dL) = planets[2]
    N = dL * T
    while (N >= 360.0):
      N = N - 360.0
    while (N < 0.0):
      N = N + 360.0
# FIXME: Should adjust the orbit's inclination and orientation
#    I = I0 + dI * T
#    wbar = wbar0 + dwbar * T
#    O = O0 + dO * T
    M = N + L0 - wbar
    while (M >= 360.0):
      M = M - 360.0
    while (M < 0.0):
      M = M + 360.0
    Le = N + 2 * rtod(e) * sin(dtor(M)) + L0
    while (Le >= 360.0):
      Le = Le - 360.0
    while (Le < 0.0):
      Le = Le + 360.0
    v = Le - wbar
    Re = a * (1 - e * e) / (1 + e * cos(dtor(v)))

# now plot the planets.

    for i in range(len(planets)):
      if (i == 2):
        continue
      (name, wbar, e, a, I, O, L0, dL) = planets[i]
      N = dL * T
      while (N >= 360.0):
        N = N - 360.0
      while (N < 0.0):
        N = N + 360.0
# FIXME:  For Mercury, especially, dI, dwbar and dO are large enough to matter.
#      I = I0 + dI * T
#      wbar = wbar0 + dwbar * T
#      O = O0 + dO * T

# compute the mean anomaly

      M = N + L0 - wbar
      while (M >= 360.0):
        M = M - 360.0
      while (M < 0.0):
        M = M + 360.0

# compute the heliocentric longitude
# FIXME: Rather than solve Keppler's equation (M = E - 180/pi * e * sin(E)) for
#        E, we will use the approximation E ~ M and calculate the heliocentric
#        longitude using the mean anomaly instead of the eccentric anomaly.
#        This approximation is pretty close except for Mercury and Pluto which
#        have values of e high enough to make a difference.

      l = N + 2 * rtod(e) * sin(dtor(M)) + L0
      while (l >= 360.0):
        l = l - 360.0
      while (l < 0.0):
        l = l + 360.0

# now calculate the actual anomaly

      v = l - wbar

# find the planet's radial distance

      r = a * (1 - e * e) / (1 + e * cos(dtor(v)))

# calculate the heliocentric latitude

      psi = rtod(asin(sin(dtor(l - O)) * sin(dtor(I))))

# project to the ecliptic

      y = sin(dtor(l - O)) * cos(dtor(I))
      x = cos(dtor(l - O))
      k = rtod(atan2(y, x))
      lprime = k + O
      while (lprime >= 360.0):
        lprime = lprime - 360.0
      while (lprime < 0.0):
        lprime = lprime + 360.0
      rprime = r * cos(dtor(psi))

# Using the coordinates we already calculated for the current position of the
# earth, convert the above coordinates to geocentric latitude and longitude

      if (i < 2):

# for inner planets, use this equation to get longitude:

        y = rprime * sin(dtor(Le - lprime))
        x = Re - rprime * cos(dtor(Le - lprime))
        k = rtod(atan2(y, x))
        lam = 180.0 + Le + k 
      else:

# for outer planets, use this equation to get longitude:

        y = Re * sin(dtor(lprime - Le))
        x = rprime - Re * cos(dtor(lprime - Le))
        k = rtod(atan2(y, x))
        lam = k + lprime
      while (lam >= 360.0):
        lam = lam - 360.0
      while (lam < 0.0):
        lam = lam + 360.0

# all planets use the same equation for ecliptic latitude (and this atan has no quadrant ambiguity):

      y = rprime * tan(dtor(psi)) * sin(dtor(lam - lprime))
      x = Re * sin(dtor(lprime - Le))
      beta = rtod(atan(y/x))

# convert lam and beta to RA and DEC

      y = sin(dtor(lam)) * cos(eps) - tan(dtor(beta)) * sin(eps)
      x = cos(dtor(lam))
      k = rtod(atan2(y, x))
      ra = k / 15.0
      while (ra < 0.0):
        ra = ra + 24.0
      dec = rtod(asin(sin(dtor(beta)) * cos(eps) + cos(dtor(beta)) * \
                      sin(eps) * sin(dtor(lam))))

# convert to azalt

      azalt = polartoazalt((ra, dec), now)

# convert to x,y and plot if the planet is above the horizon

      if (azalt[1] >= 0.0175):
        (px, py) = self.azalttoxy(azalt)
        px = px + self.margin - 2 + self.xoffset
        py = py + self.margin - 2 + self.yoffset
        self.plot_planetary_symbol(i, px, py)

# Plot the sun.  This is virtually the same as for a planet, but we can simplify
# some computations because the sun is by definition on the ecliptic.

    (name, wbar, e, a, I, O, L0, dL) = sun
    N = dL * T
    while (N >= 360.0):
      N = N - 360.0
    while (N < 0.0):
      N = N + 360.0    
    M = N + L0 - wbar
    while (M >= 360.0):
      M = M - 360.0
    while (M < 0.0):
      M = M + 360.0
    v = M + 2 * rtod(e) * sin(dtor(M))
    while (v >= 360.0):
      v = v - 360.0
    while (lam < 0.0):
      v = v + 360.0
    lam = v + wbar
    while (lam >= 360.0):
      lam = lam - 360.0
    while (lam < 0.0):
      lam = lam + 360.0
    y = sin(dtor(lam)) * cos(eps)
    x = cos(dtor(lam))
    k = rtod(atan2(y, x))
    ra = k / 15.0
    while (ra < 0.0):
      ra = ra + 24.0

# because beta is (by definition) 0, calculating dec is far simpler:

    dec = rtod(asin(sin(sin(eps) * sin(dtor(lam)))))
    azalt = polartoazalt((ra, dec), now)
    if (azalt[1] >= 0.0175):
      (px, py) = self.azalttoxy(azalt)
      px = px + self.margin - 2 + self.xoffset
      py = py + self.margin - 2 + self.yoffset
      self.plot_planetary_symbol(7, px, py)

# Plot the moon.  Since the moon's orbit varies radically over time, there are
# a lot of "fudge factor" correction terms in these computations.

    (name, L0, P0, N0, I, e, a, phi0, tau) = moon
    D = jtime(now) - float(tau)
    N = 360.0 / 365.242191 * D
    while (N >= 360.0):
      N = N - 360.0
    while (N < 0.0):
      N = N + 360.0
    M = N + L0 - P0
    while (M >= 360.0):
      M = M - 360.0
    while (M < 0.0):
      M = M + 360.0
    Ec = 2 * rtod(e) * sin(dtor(M))
    lam = Ec + L0
    while (lam >= 360.0):
      lam = lam - 360.0
    while (lam < 0.0):
      lam = lam + 360.0
    l = 13.176396 * D + L0               
    while (l >= 360.0):
      l = l - 360.0
    while (l < 0.0):
      l = l + 360.0
    Mm = l - 0.1114041 * D - P0
    while (Mm >= 360.0):
      Mm = Mm - 360.0
    while (Mm < 0.0):
      Mm = Mm + 360.0
    N = N0 - 0.0529539 * D
    while (N >= 360.0):
      N = N - 360.0
    while (N < 0.0):
      N = N + 360.0
    Ev = 1.2739 * sin(dtor(2 * (1 - lam) - Mm))
    Ae = 0.1858 * sin(dtor(M))
    A3 = 0.37 * sin(dtor(M))
    Mprime = Mm + Ev - Ae - A3
    Ec = 6.2886 * sin(dtor(Mprime))
    A4 = 0.214 * sin (dtor(2 * Mprime))
    lprime = l + Ev  + Ec - Ae + A4
    V = 0.6583 * sin(dtor(2 * (lprime - lam)))
    lon = lprime + V
    Nprime = N - 0.16 * sin(dtor(M))
    y = sin(dtor(lon - Nprime)) * cos(dtor(I))
    x = cos(dtor(lon - Nprime))
    k = rtod(atan2(y, x))
    lamM = Nprime + k
    while (lamM >= 360.0):
      lamM = lamM - 360.0
    while (lamM < 0.0):
      lamM = lamM + 360.0
    betM = rtod(asin(sin(dtor(lon - Nprime)) * sin(dtor(I))))
    y = sin(dtor(lamM)) * cos(eps) - tan(dtor(betM)) * sin(eps)
    x = cos(dtor(lamM))
    k = rtod(atan2(y, x))
    ra = k / 15.0
    while (ra < 0.0):
      ra = ra + 24.0
    while (ra > 24.0):
      ra = ra - 24.0
    dec = rtod(asin(sin(dtor(betM)) * cos(eps) + cos(dtor(betM)) * \
                    sin(eps) * sin(dtor(lamM))))
# FIXME: We aren't accounting for parallax: RA and DEC are actual, not apparent.
#	This could introduce as much as a degree (approx.) of error.

# convert to azalt

    azalt = polartoazalt((ra, dec), now)
    if (azalt[1] >= 0.0175):
      (px, py) = self.azalttoxy(azalt)
      px = px + self.margin - 2 + self.xoffset
      py = py + self.margin - 2 + self.yoffset
      self.plot_planetary_symbol(2, px, py)
    self.gc.set_foreground(self.colors[1])
    return True


  def callback(self, widget, data=None):

# Handle control changes here.

    global nightvision
    global invertdisplay
    global fliphorizontally
    global drawconstellations
    global limitingmagnitude
    global longitude
    global latitude
    global specifytime
    global lon_alert
    global lat_alert
    global time_alert
    global zone_alert

    if (data == "night vision"):
      nightvision = button1.get_active()
      self.plotchart()
    if (data == "invert display"):
      invertdisplay = button2.get_active()
      self.plotchart()
    if (data == "flip horizontally"):
      fliphorizontally = button3.get_active()
      self.plotchart()
    if (data == "draw constellations"):
      drawconstellations = button4.get_active()
      self.plotchart()
    if (data == "location change"):
      s = entry1.get_text()
      if (s.find("m") >= 0):
        lon = angletofloat(s)
      else:
        try:
          lon = float(s)
        except:
          lon = -1.0
      if ((lon < 0.0) or (lon > 180.0)):
# FIXME: beep
	lon_alert.show()
        entry1.set_text(floattoangle(abs(longitude)))
        if (longitude < 0.0):
          rb2.clicked()
        else:
          rb1.clicked()
        return True
      else:
        lon_alert.hide()
        longitude = lon
        entry1.set_text(floattoangle(abs(longitude)))
        if (rb2.get_active()):
          longitude = -longitude
      s = entry2.get_text()
      if (s.find("m") >= 0):
        lat = angletofloat(s)
      else:
        try:
          lat = float(s)
        except:
          lat = -1.0
      if ((lat < 0.0) or (lat > 90.0)):
# FIXME: beep
        lat_alert.show()
        entry2.set_text(floattoangle(abs(latitude)))
        if (latitude < 0.0):
          rb4.clicked()
        else:
          rb3.clicked()
        return True
      else:
        lat_alert.hide()
        latitude = lat
        entry2.set_text(floattoangle(abs(latitude)))
        if (rb4.get_active()):
          latitude = -latitude
      self.plotchart()
# FIXME: error reporting for time and zone does not appear to work but it will automatically reset the values to "now" and "local_zone" in event of error.
    if (data == _("time change")):
      specifytime = rb6.get_active()
      if (specifytime):
        if (syntax_check_time() == False):
          time_alert.show()
          specifytime = False
          now = datetime.utcnow()
          (tstr, ostr) = set_time_and_UTC_offset()
          entry3.set_text(tstr)
          entry4.set_text(ostr)
          return True
        else:
          time_alert.hide()
          if (syntax_check_zone() == False):
            zone_alert.show()
            specifytime = False
            now = datetime.utcnow()
            (tstr, ostr) = set_time_and_UTC_offset()
            entry3.set_text(tstr)
            entry4.set_text(ostr)
            return True
          else:
            zone_alert.hide()
      self.plotchart()
    if (data == _("rb7 clicked")):
      limitingmagnitude = 0.0
      self.plotchart()
    if (data == _("rb8 clicked")):
      limitingmagnitude = 1.0
      self.plotchart()
    if (data == _("rb9 clicked")):
      limitingmagnitude = 2.0
      self.plotchart()
    if (data == _("rb10 clicked")):
      limitingmagnitude = 3.0
      self.plotchart()
    if (data == _("rb11 clicked")):
      limitingmagnitude = 4.0
      self.plotchart()
    if (data == _("rb12 clicked")):
      limitingmagnitude = 5.0
      self.plotchart()
    if (data == _("rb13 clicked")):
      limitingmagnitude = 6.0
      self.plotchart()

      
  def cleararea(self):
    
# Clear the drawing surface

    if (nightvision):
      self.gc.set_foreground(self.colors[1])
    else:
      self.gc.set_foreground(self.colors[3])
    self.window.draw_rectangle(self.gc,
                                    True,
                                    1,
                                    1,
                                    self.screensize[0],
                                    self.screensize[1]);
    label1.queue_draw()
    label2.queue_draw()
    label4.queue_draw()
    label5.queue_draw()
    label6.queue_draw()
    button1.queue_draw()
    button2.queue_draw()
    button3.queue_draw()
    button4.queue_draw()
    button5.queue_draw()
    button6.queue_draw()
    rb1.queue_draw()
    rb2.queue_draw()
    rb3.queue_draw()
    rb4.queue_draw()
    rb5.queue_draw()
    rb6.queue_draw()
    rb7.queue_draw()
    rb8.queue_draw()
    rb9.queue_draw()
    rb10.queue_draw()
    rb11.queue_draw()
    rb12.queue_draw()
    rb13.queue_draw()


  def plot_planetary_symbol(self, i, px, py):

# i is planet number (0 = mercury, 1= venus, 2 = moon, 3 = mars, 4 = jupiter,
# 5 = saturn, 6 = uranus, 7 = sun); (px, py) is the center point of the symbol.

    if (i == 0):

# mercury

      self.gc.set_line_attributes(2, gtk.gdk.LINE_SOLID, gtk.gdk.CAP_BUTT,
                                  gtk.gdk.JOIN_MITER)
      self.window.draw_arc(self.gc, False, px-12, py-12, 24, 24, 0, 360*64)
      self.window.draw_arc(self.gc, False, px-5, py-7, 10, 10, 0, 360*64)
      self.window.draw_line(self.gc, px+4, py-9, px+4, py-7)
      self.window.draw_line(self.gc, px-4, py-9, px-4, py-7)
      self.gc.set_line_attributes(1, gtk.gdk.LINE_SOLID, gtk.gdk.CAP_BUTT,
                                  gtk.gdk.JOIN_MITER)
      self.window.draw_line(self.gc, px, py+3, px, py+7)
      self.window.draw_line(self.gc, px-2, py+5, px+2, py+5)
    elif (i == 1):

# venus

      self.gc.set_line_attributes(2, gtk.gdk.LINE_SOLID, gtk.gdk.CAP_BUTT,
                                  gtk.gdk.JOIN_MITER)
      self.window.draw_arc(self.gc, False, px-12, py-12, 24, 24, 0, 360*64)
      self.window.draw_arc(self.gc, False, px-5, py-7, 10, 10, 0, 360*64)
      self.gc.set_line_attributes(1, gtk.gdk.LINE_SOLID, gtk.gdk.CAP_BUTT,
                                  gtk.gdk.JOIN_MITER)
      self.window.draw_line(self.gc, px, py+3, px, py+7)
      self.window.draw_line(self.gc, px-2, py+5, px+2, py+5)
    elif (i == 2):

# moon

      self.gc.set_line_attributes(2, gtk.gdk.LINE_SOLID, gtk.gdk.CAP_BUTT,
                                  gtk.gdk.JOIN_MITER)
      self.window.draw_arc(self.gc, False, px-12, py-12, 24, 24, 0, 360*64)
      self.window.draw_polygon(self.gc, True, ((px+1, py-11), (px+4, py-11),
                                               (px+5, py-10), (px+6, py-9),
                                               (px+7, py-8),  (px+8, py-7),
                                               (px+10, py-2), (px+12,py),
                                               (px+10, py+2), (px+8, py+7),
                                               (px+7, py+8),  (px+6, py+9),
                                               (px+5, py+10), (px+4,py+11),
                                               (px+1, py+11), (px+4, py+4),
                                               (px+6, py+2),  (px+6, py-2),
                                               (px+4, py-4)))
      self.gc.set_line_attributes(1, gtk.gdk.LINE_SOLID, gtk.gdk.CAP_BUTT,
                                  gtk.gdk.JOIN_MITER)
    elif (i == 3):

# mars

      self.gc.set_line_attributes(2, gtk.gdk.LINE_SOLID, gtk.gdk.CAP_BUTT,
                                  gtk.gdk.JOIN_MITER)
      self.window.draw_arc(self.gc, False, px-12, py-12, 24, 24, 0, 360*64)
      self.window.draw_arc(self.gc, False, px-6, py-4, 10, 10, 0, 360*64)
      self.gc.set_line_attributes(1, gtk.gdk.LINE_SOLID, gtk.gdk.CAP_BUTT,
                                  gtk.gdk.JOIN_MITER)
      self.window.draw_line(self.gc, px+2, py-2, px+6, py-6)
      self.window.draw_line(self.gc, px+3, py-6, px+6, py-6)
      self.window.draw_line(self.gc, px+6, py-6, px+6, py-3)
    elif (i == 4):

# jupiter

      self.gc.set_line_attributes(2, gtk.gdk.LINE_SOLID, gtk.gdk.CAP_BUTT,
                                  gtk.gdk.JOIN_MITER)
      self.window.draw_arc(self.gc, False, px-12, py-12, 24, 24, 0, 360*64)
      self.window.draw_line(self.gc, px-6, py-6, px-4, py-8)
      self.window.draw_line(self.gc, px-4, py-8, px-2, py-8)
      self.window.draw_line(self.gc, px-2, py-8, px+1, py-6)
      self.window.draw_line(self.gc, px+1, py-6, px-5, py+2)
      self.window.draw_line(self.gc, px-5, py+2, px+7, py+2)
      self.window.draw_line(self.gc, px+4, py-8, px+4, py+7)
      self.gc.set_line_attributes(1, gtk.gdk.LINE_SOLID, gtk.gdk.CAP_BUTT,
                                  gtk.gdk.JOIN_MITER)
    elif (i == 5):

# saturn

      self.gc.set_line_attributes(2, gtk.gdk.LINE_SOLID, gtk.gdk.CAP_BUTT,
                                  gtk.gdk.JOIN_MITER)
      self.window.draw_arc(self.gc, False, px-12, py-12, 24, 24, 0, 360*64)
      self.window.draw_line(self.gc, px-6, py-6, px-6, py+5)
      self.window.draw_line(self.gc, px-6, py, px-5, py-1)
      self.window.draw_line(self.gc, px-5, py-1, px-4, py-2)
      self.window.draw_line(self.gc, px-4, py-2, px-1, py-3)
      self.window.draw_line(self.gc, px-1, py-3, px, py-4)
      self.window.draw_line(self.gc, px, py-4, px+1, py+1)
      self.window.draw_line(self.gc, px+1, py+1, px-1, py+4)
      self.window.draw_line(self.gc, px-1, py+4, px, py+5)
      self.window.draw_line(self.gc, px, py+5, px+6, py+4)
      self.gc.set_line_attributes(1, gtk.gdk.LINE_SOLID, gtk.gdk.CAP_BUTT,
                                  gtk.gdk.JOIN_MITER)
    elif (i == 6):

# uranus

      self.gc.set_line_attributes(2, gtk.gdk.LINE_SOLID, gtk.gdk.CAP_BUTT,
                                  gtk.gdk.JOIN_MITER)
      self.window.draw_arc(self.gc, False, px-12, py-12, 24, 24, 0, 360*64)
      self.window.draw_arc(self.gc, False, px-5, py-3, 10, 10, 0, 360*64)
      self.window.draw_arc(self.gc, True, px-2, py, 4, 4, 0, 360*64)
      self.gc.set_line_attributes(1, gtk.gdk.LINE_SOLID, gtk.gdk.CAP_BUTT,
                                  gtk.gdk.JOIN_MITER)
      self.window.draw_line(self.gc, px, py-3, px, py-9)
      self.window.draw_line(self.gc, px-2, py-5, px, py-9)
      self.window.draw_line(self.gc, px+2, py-5, px, py-9)
    else:

# sun

      self.gc.set_line_attributes(2, gtk.gdk.LINE_SOLID, gtk.gdk.CAP_BUTT,
                                  gtk.gdk.JOIN_MITER)
      self.window.draw_arc(self.gc, False, px-12, py-12, 24, 24, 0, 360*64)
      self.window.draw_arc(self.gc, True, px-2, py-2, 4, 4, 0, 360*64)
      self.gc.set_line_attributes(1, gtk.gdk.LINE_SOLID, gtk.gdk.CAP_BUTT,
                                  gtk.gdk.JOIN_MITER)


  def plot_DSO(self, type, maja, mina, mag, px, py):
    if (not invertdisplay):
      if (nightvision):
        fg_color = self.colors[2]
      else:
        fg_color = self.colors[0]
    else:
      fg_color = self.colors[1]
    if (mag > 9.0) or (maja < 10.0):
# this object is too dim or too small to be interesting
      pass;
    else:
# FIXME:  The Messier Catalog should include the object's position angle and this
#         code should plot it in the proper orientation.  What we have here will
#         always plot the object with the major axis left-right and the minor
#         axis up-down.

# We plot these objects scaled up by a decreasing amount but always much larger
# than actual size.  Otherwise, they would be mostly too small to plot.
      if (maja < 20.0):
        maja = maja * 8.0
        mina = mina * 8.0
      elif (maja < 40.0):
        maja = maja * 6.0
        mina = mina * 6.0
      elif (maja < 80.0):
        maja = maja * 5.0
        mina = mina * 5.0
      elif (maja < 120.0):
        maja = maja * 4.0
        mina = mina * 4.0
      else:
        maja = maja * 3.0
        mina = mina * 3.0
      dx = maja / 60.0 * self.diameter / 180.0
      dy = mina / 60.0 * self.diameter / 180.0
      if (type == 'Gal'):
# plot as gray ellipse with solid outline.
        self.gc.set_foreground(self.colors[3])
        self.window.draw_arc(self.gc,
                              True,
                              px - dx / 2,
                              py - dy / 2,
                              dx,
                              dy,
                              0,
                              360 * 64)
        self.gc.set_foreground(fg_color)
        self.window.draw_arc(self.gc,
                              False,
                              px - dx / 2,
                              py - dy / 2,
                              dx,
                              dy,
                              0,
                              360 * 64)
      elif (type == 'PlN'):
# plot as gray circle with central dot
        self.gc.set_foreground(self.colors[3])
        self.window.draw_arc(self.gc,
                              True,
                              px - dx / 2,
                              py - dx / 2,
                              dx,
                              dx,
                              0,
                              360 * 64)
        self.gc.set_foreground(fg_color)
        self.window.draw_arc(self.gc,
                              True,
                              px - 2,
                              py - 2,
                              4,
                              4,
                              0,
                              360 * 64)
      elif (type == 'SNR') or (type == 'OCl'):
# plot as gray circle with no outline.
        self.gc.set_foreground(self.colors[3])
        self.window.draw_arc(self.gc,
                              True,
                              px - dx / 2,
                              py - dx / 2,
                              dx,
                              dx,
                              0,
                              360 * 64)
        self.gc.set_foreground(fg_color)
      elif (type == 'C/N') or (type == 'DfN'):
# plot as gray rectangle with no outline.
        self.gc.set_foreground(self.colors[3])
        self.window.draw_rectangle(self.gc,
                                        True,
                                        px - dx / 2,
                                        py - dy / 2,
                                        dx,
                                        dy);
        self.gc.set_foreground(fg_color)
      elif (type == 'GCl'):
# plot as gray circle with outline and central dot.
        self.gc.set_foreground(self.colors[3])
        self.window.draw_arc(self.gc,
                              True,
                              px - dx / 2,
                              py - dx / 2,
                              dx,
                              dx,
                              0,
                              360 * 64)
        self.gc.set_foreground(fg_color)
        self.window.draw_arc(self.gc,
                              False,
                              px - dx / 2,
                              py - dx / 2,
                              dx,
                              dx,
                              0,
                              360 * 64)
        self.window.draw_arc(self.gc,
                              True,
                              px - 2,
                              py - 2,
                              4,
                              4,
                              0,
                              360 * 64)
      else:
#	Dbl = double star
#	??? = unknown or unclassified object
# these are not plotted.
        pass

# ============================== StarChart Object ===============================

class StarChart(activity.Activity):
  def __init__(self, handle):
    global now
    global zoneoffset
    activity.Activity.__init__(self, handle)
    os.chdir(get_bundle_path())
    self.set_title(_("Star Chart Activity"))
                    
# Iniitialize time to now and offset to our zone.

    now = datetime.utcnow()
    (tstr, ostr) = set_time_and_UTC_offset()
    (hh, mm) = parse_zone_offset(ostr)
    zoneoffset = 60 * hh
    if (hh < 0):
      zoneoffset = zoneoffset - mm
    else:
      zoneoffset = zoneoffset + mm

# Create toolbox
      
    toolbox = activity.ActivityToolbox(self)
    self.set_toolbox(toolbox)

    self.what_toolbar = gtk.Toolbar();
    self.where_toolbar = gtk.Toolbar();
    self.when_toolbar = gtk.Toolbar();
      
# Fill the toolbox bars

    self.what_toolbar.add(button1)
    self.what_toolbar.add(button2)
    self.what_toolbar.add(button3)
    self.what_toolbar.add(button4)
    self.what_toolbar.add(label6)
    self.what_toolbar.add(rb7)
    self.what_toolbar.add(rb8)
    self.what_toolbar.add(rb9)
    self.what_toolbar.add(rb10)
    self.what_toolbar.add(rb11)
    self.what_toolbar.add(rb12)
    self.what_toolbar.add(rb13)
    self.where_toolbar.add(label1)
    self.where_toolbar.add(entry1)
    self.where_toolbar.add(rb1)
    self.where_toolbar.add(rb2)
    self.where_toolbar.add(label2)
    self.where_toolbar.add(entry2)
    self.where_toolbar.add(rb3)
    self.where_toolbar.add(rb4)
    self.where_toolbar.add(button5)
    self.when_toolbar.add(rb5)
    self.when_toolbar.add(rb6)
    self.when_toolbar.add(label4)
    self.when_toolbar.add(entry3)
    self.when_toolbar.add(label5)
    self.when_toolbar.add(entry4)
    self.when_toolbar.add(button6)
    toolbox.add_toolbar('What', self.what_toolbar)
    toolbox.add_toolbar('Where', self.where_toolbar)
    toolbox.add_toolbar('When', self.when_toolbar)
    scrolled = gtk.ScrolledWindow()
    scrolled.set_policy(gtk.POLICY_NEVER, gtk.POLICY_AUTOMATIC)
    scrolled.props.shadow_type = gtk.SHADOW_NONE
    self.chart = ChartDisplay()
    eb = gtk.EventBox()
    eb.add(self.chart)
    eb.modify_bg(gtk.STATE_NORMAL, gtk.gdk.color_parse("gray"))
    scrolled.add_with_viewport(eb)

# create the alerts

    lon_alert.props.title = _("Input Error")
    lon_alert.props.msg = _("The longitude was not understood.")
    self.add_alert(lon_alert)
    lat_alert.props.title = _("Input Error")
    lat_alert.props.msg = _("The latitude  was not understood.")
    self.add_alert(lat_alert)
    time_alert.props.title = _("Input Error")
    time_alert.props.msg = _("The time was not understood.")
    self.add_alert(time_alert)
    zone_alert.props.title = _("Input Error")
    zone_alert.props.msg = _("The time zone offset was not understood.")
    self.add_alert(zone_alert)


# Connect the events

    button1.connect("clicked", self.chart.callback, "night vision")
    button2.connect("clicked", self.chart.callback, "invert display")
    button3.connect("clicked", self.chart.callback, "flip horizontally")
    button4.connect("clicked", self.chart.callback, "draw constellations")
    rb7.connect("clicked", self.chart.callback, "rb7 clicked")
    rb8.connect("clicked", self.chart.callback, "rb8 clicked")
    rb9.connect("clicked", self.chart.callback, "rb9 clicked")
    rb10.connect("clicked", self.chart.callback, "rb10 clicked")
    rb11.connect("clicked", self.chart.callback, "rb11 clicked")
    rb12.connect("clicked", self.chart.callback, "rb12 clicked")
    rb13.connect("clicked", self.chart.callback, "rb13 clicked")
    button5.connect("clicked", self.chart.callback, "location change")
    button6.connect("clicked", self.chart.callback, "time change")
    self.chart.connect("expose_event", self.chart.area_expose_cb)
    toolbox.show()
    self.chart.show()
    eb.show()
    scrolled.show()
    self.set_canvas(scrolled)
    scrolled.show()
    self.show_all()
    toolbar = toolbox.get_activity_toolbar()
    toolbar.share.hide()
    initialize_controls()
    self.chart.plotchart()
    lon_alert.hide()

  def read_file(self, filename):
    global nightvision
    global invertdisplay
    global fliphorizontally
    global drawconstellations
    global limitingmagnitude
    global latitude
    global longitude
    global specifytime
    global zoneoffset
    global now

    f = open(filename, "r")
    nightvision = bool(int(self.metadata.get('Night_Vision', '0')))
    invertdisplay = bool(int(self.metadata.get('Invert', '0')))
    fliphorizontally = bool(int(self.metadata.get('Flip', '0')))
    drawconstellations = bool(int(self.metadata.get('Constellations', '1')))
    limitingmagnitude = float(self.metadata.get('Magnitude', '4.0'))
    latitude =  float(self.metadata.get('Latitude', '42.64333333'))
    longitude = float(self.metadata.get('Longitude', '-71.3963888889'))
    specifytime = bool(int(self.metadata.get('Specify_Time', '0')))
    if (specifytime):
      ts = self.metadata.get('Time', now.strftime("%Y/%m/%d,%H:%M"))
      zs = self.metadata.get('Zone_Offset', str(zoneoffset))
      entry3.set_text(ts)
      entry4.set_text(zs)
      now = get_time_and_UTC_offset(entry3.get_text(), entry4.get_text())
      (hh, mm) = parse_zone_offset(entry4.get_text())
      zoneoffset = 60 * hh
      if (hh < 0):
        zoneoffset = zoneoffset - mm
      else:
        zoneoffset = zoneoffset + mm
    f.close()
    initialize_controls()
    self.chart.plotchart()
    
    
  def write_file(self, filename):
    f = open(filename, "w")
    self.metadata['Night_Vision'] = str(int(nightvision))
    self.metadata['Invert'] = str(int(invertdisplay))
    self.metadata['Flip'] = str(int(fliphorizontally))
    self.metadata['Constellations'] = str(int(drawconstellations))
    self.metadata['Magnitude'] = str(limitingmagnitude)
    self.metadata['Latitude'] = str(latitude)
    self.metadata['Longitude'] = str(longitude)
    self.metadata['Specify_Time'] = str(int(specifytime))
# Unlike the other settings, it's easier to store the time and zone offset as
# they are represented in the text entry controls than to attempt to convert
# now and zoneoffset into a representation of local time and offset.
    self.metadata['Time'] = entry3.get_text()
    self.metadata['Zone_Offset'] = entry4.get_text()
    f.close()

