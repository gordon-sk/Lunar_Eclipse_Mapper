# Midshipman 1/C Gordon Kiesling, Vanderbilt University
# Naval Research Enterprise Internship Program summer undergraduate intern
# United States Naval Observatory, Washington DC
# Under direction of Dr. Malynda Frouard, USNO
#
# This program emulates the lunar eclipse plotting as done by
# Her Majesty's Nautical Almanac Office in Fortran
# It requires user input of a year and ∆T value and will present them with
# the eclipses that occur that year and maps for the specified lunar eclipses
# All effort has been made to emulate the fortran plots in aesthetic as well
# as data accuracy.
#
# This program was written in Python 3.6.1 as installed by Anaconda.
# It requires the libraries specified below, most of which are either native
# to Python or come with anaconda.
# Basemap may require an install.
#
# To run, this program also requires three images to be placed in the same
# directory as this file -- compass.png, big_moon.png, and big_moon_red.png
# These file names are declared at the top of the file, directly below the
# imports, for ease of modification if necessary.
#
# It also requires high-accuracy ephemeris data. The program was written
# with 405 data downloaded from NASA's public file directory located at
# ftp://ssd.jpl.nasa.gov/pub/eph/planets/ascii/de405/
# The ephemeris file is also named at the top of the program
#
#
#
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from matplotlib import text as mtext
from matplotlib.offsetbox import (OffsetImage, AnnotationBbox)
from matplotlib.cbook import get_sample_data
from PIL import Image
import novas.compat as nov
import novas.compat.eph_manager as eph_mgr
import novas.constants as const
from math import radians as rad
from math import sin
from math import cos
from math import fabs
import numpy as np
import math
import calendar
from sys import path

# time units are Julian days or fractions thereof
initial_sensitivity = 2
initial_time_step = 2 * pow(10, -initial_sensitivity)    # 28.8 minutes, in days
map_t_step = 1 / 24 / 60           #  1 minute

# Inclination of the moon's orbit to the eccliptic, in degrees
I = 5.1453964
# Corrections for difference between moon's visible center and center of mass
# both in units of degrees
α_correction = .5 / 60 / 60
δ_correction = -.25 / 60 / 60
# Celestial obejects -- for use by novas
moon = nov.make_object(0, 11, "moon", None)
sun = nov.make_object(0, 10, "sun", None)
# Ephemeris data file name
Ephem_file = "JPLEPH405"
# Image file names
compass_img_file = path[0] + "/compass.png"
moon_img_file = path[0] + "/big_moon.png"
red_moon_img_file = path[0] + "/big_moon_red.png"


# Main method. Calls the helper methods, acts as table of contents.
# Prints useful information to the user.
# -------------
# Input: none
# Output: none
#
def main():
    print(" This program lets the user specify a year in which eclipses take\n",
      "place, and lets them choose specifics eclipses that are available.\n",
      "It then outputs data about the said eclipse, including maps, dates,\n",
      "& times. All data is pulled from NASA's publicly available ephemeris\n",
      "files located at ftp://ssd.jpl.nasa.gov/pub/eph/planets/ascii/",
      "\n Specific Ephemeris Data: LE405/DE405\n")
    print(" All dates in (Y, M, D, time) format UTC\n")
    ΔT = dT_prompt()                            # User chooses ∆t value to use
    dates = eph_mgr.ephem_open(Ephem_file)      # for printing range of data
    year = year_prompt(dates, )                   # User selects year
    eclipses = eclipse_finder(year)             # Program searches for eclipses
    which_event = event_prompt(eclipses)        # User chooses eclipse to map

    # a for loop is the tighest way to write this, considering
    # the 'all' option for mapping
    for n in range(eclipses.__len__()):
        lunar = eclipses[n][2]
        if which_event - 1 == n or which_event == 'all':
            if lunar:
                lunar_eclipse_mapper(eclipses, n, ΔT)
            else:
                solar_eclipse_mapper(eclipses, n, ΔT)

    # Close up the memory
    eph_mgr.ephem_close()
    # end main


# Method year_prompt
# Prompts the user for a year to search for eclipses in.
# Checks date to ensure it is both a valid object type and is within the range
# of data provided by the ephemeris file.
# Continues to prompt user until date is correct
# -------------
# Input: User types year of interest
# Output: Year is returned to main for use
#
def year_prompt(dates):
    print(" Enter the year in which the eclipse you are" +
                 " interested in is occuring (Gregorian calendar): ")
    year = input(" Range of data is " + str(nov.cal_date(dates[0])) + " to " +
                 str(nov.cal_date(dates[1])) + ": ")
    try:
        year = int(year)
        if year < nov.cal_date(dates[0])[0] or year > nov.cal_date(dates[1])[0]:
            print("\n Year is out of range. Try again.\n\n")
            year = year_prompt(dates)
    except ValueError:
        print("\n Input is not an integer. Try again.\n\n")
        year = year_prompt(dates)
    return year


# Method ∆T_prompt
# Prompts user for a d∆ value for correcting between TT time and UTC time.
# Checks to ensure input is a number. No checks to ensure number makes sense.
# -------------
# Input: User inputs either integer or β_mimal number for ∆T value
# Output: ∆T is converted from seconds to Julian days and returned to main
#
def dT_prompt():
    ΔT = input(" Enter the delta T value to use (seconds): ")
    try:
        ΔT = float(ΔT)
    except ValueError:
        print("\n Input is not a number. Try again.\n\n")
        ΔT = dT_prompt()
    print()
    return ΔT / (60 * 60 * 24)


# Method event_prompt
# Prints eclipses found by method eclipse_mapper in chronological order,
# and has user select eclipse of interest
# User can also choose to map 'all' eclipses. This was implemented to save
# some time for anyone doing a lot of mapping.
# -------------
# Input: user specifies integer to indicate which eclipse to map, or chooses
#  to map all 'all' input is controlled for extra spaces, capital letters
# Output: variable which_event returned for use by mapper
#
def event_prompt(eclipses):
    print("\n  ", eclipses.__len__(), "eclipses found:")
    print("   Note: Hours are aproximate")
    for x in range(eclipses.__len__()):
        date = nov.cal_date(eclipses[x][1])
        date = (date[0], date[1], date[2], round(date[3], 1))
        print("  ", x + 1, ": ", eclipses[x][0], "occuring on", date)
    print("\n Please enter the number corresponding to the eclipse of interest.")
    which_event = input(" Alternatively, type 'all' to receive output" +
                        " for all listed events: ")
    try:
        which_event = int(which_event)
        if which_event > eclipses.__len__() or which_event < 1:
            print("\n Input is out of range. Try again.\n\n")
            which_event = event_prompt(eclipses)
    except ValueError:
        if which_event.lower() == "all":
            which_event = which_event.lower().replace(" ", "")
        else:
            print("\n Input is invalid. Try again.\n\n")
            which_event = event_prompt(eclipses)
    return which_event


# method eclipse_finder
# Does all the math and time iterating to find eclipses.
# Timestep is specified before main, below imports, and default is 28.8 minutes
# or .02 Julian days. All math is taken from the explanatory supplement,
# pages 457-464.
# -------------
# Input: integer year is passed from main
# Output: A list of eclipses, along with the Julian date at which they occur
# and a boolean indicating whether or not they are lunar (True) or solar (False)
#
def eclipse_finder(year):
    # dates to begin and end search
    date = nov.julian_date(year, 1, 1)
    end = nov.julian_date(year + 1, 1, 1)

    # empty list to which found eclipses will be appended
    eclipses = []

    # lists to save specific data points chronologically for later use
    angular_separations = []
    apparent_semidiameters_sun, apparent_semidiameters_moon = [], []
    parallaxes_sun, parallaxes_moon = [], []

    print(" Searching specified year for eclipses...")
    while date <= end:

        # using novas to get positional data.
        # Converting RA from hour angle to degrees.
        pos_moon = nov.app_planet(date, moon)
        pos_sun = nov.app_planet(date, sun)
        pos_moon = (pos_moon[0] * 15.0, pos_moon[1], pos_moon[2])
        pos_sun = (pos_sun[0] * 15.0, pos_sun[1], pos_sun[2])

        # Potential eclipses are only tested for if the bodies are close
        # to conjuncton or opposition. 1.5 degrees is chosen because it is
        # roughly the size of the earth's penumbra, and then some, and well
        # above the combined average radii of the sun and moon.
        alignment = ""
        if fabs(pos_moon[0] - pos_sun[0]) < 1.5:
            alignment = "conjunction"
            difference = fabs(pos_moon[0] - pos_sun[0])
        elif fabs(pos_moon[0] - pos_sun[0]) - 180 < 1.5:
            alignment = "opposition"
            difference = fabs(pos_moon[0] - pos_sun[0]) - 180
#        elif fabs(pos_sun[0] - pos_moon[0] - 180) < 1.5:
#            alignment = "opposition"
#            difference = fabs(pos_sun[0] - pos_moon[0] - 180)

        # This is all skipped if the bodies are in different parts of the sky.
        if alignment is not "":
            # Gathering variables to execute eq. 11.9 pg 459 explanatory
            # supplement, the angular separation of the sun and moon at
            # opposition or conjunction.
            # Required variables:
            # SP  = moon's motion in longitude (RA)
            # SS' = sun's motion in longitude (RA)
            # βm  = relative β_mlination of sun and moon
            # I   = angle between moon's orbit and the ecliptic
            # I is β_mlared at the top of the file, and is a constant.
            ##### sun and moon speeds in longitude ####
            # No need to include the timestep or convert to degrees,
            # as it would just divide out
            old_pos_moon = nov.app_planet(date - .00001, moon)[0]
            current_pos_moon = nov.app_planet(date, moon)[0]
            moon_vel = (current_pos_moon - old_pos_moon)

            old_pos_sun = nov.app_planet(date - .00001, sun)[0]
            current_pos_sun = nov.app_planet(date, sun)[0]
            sun_vel = (current_pos_sun - old_pos_sun)
            ######
            if alignment == "conjunction":
                βm = (pos_moon[1] - pos_sun[1])
            elif alignment == "opposition":
                βm = (pos_moon[1] + pos_sun[1])

            # Passing the variables to a helper method, which will do the math.
            σ = geocentric_ang_sep_finder(moon_vel, sun_vel, βm, date, alignment)
            # Because this equation is for bodies at opposition or conjunction,
            # and at this point there is still some longitudal separation,
            # the true angular separation is not given by eq. 11.9.
            # It is actually the hypotenuse of the triangle defined by sigma
            # and the RA diffference.
            # We redefine sigma to include this difference here
            σ[0] = math.sqrt(σ[0] ** 2 + difference ** 2)
            angular_separations.append(σ)
            ##########

            # More helper method calls for data that will be required below.
            sun_angular_semidiameter = angular_semidiameter(pos_sun[2], sun)
            moon_angular_semidiameter = angular_semidiameter(pos_moon[2], moon)
            sun_HP = parallax(pos_sun[2])
            moon_HP = parallax(pos_moon[2])

            apparent_semidiameters_sun.append(sun_angular_semidiameter)
            apparent_semidiameters_moon.append(moon_angular_semidiameter)
            parallaxes_sun.append(sun_HP)
            parallaxes_moon.append(moon_HP)

            # reset for the next iteration
            alignment = ""

        # Keeping timesteps consistent
        date = round(date + initial_time_step, initial_sensitivity)

    # Now that we have data corresponding to the bodies at conjunction or
    # opposition, we test it using the formulas 11.13, 11.14, 11.15, 11.21
    # and 11.23
    # this is not done in the loop above because we want to specify that sigma
    # is less than both the value of sigma before and after it, to find the
    # local mininums that are most likely to be eclipses
    # The type of eclipse, date of occurance, and a whether or not it is lunar
    # (for ease of sorting in main) are appended to our master eclipse list
    for x in range(angular_separations.__len__()):
        # finding the local minimum sigma in part of our dataset
        less = False
        sigma = fabs(angular_separations[x][0])
        try:
            if sigma < fabs(angular_separations[x-1][0]):
                if sigma < fabs(angular_separations[x+1][0]):
                    less = True
        except IndexError:
            pass    # catching the first and last value in the list
        if less:
            π_s = parallaxes_sun[x]
            π_m = parallaxes_moon[x]
            s_s = apparent_semidiameters_sun[x]
            s_m = apparent_semidiameters_moon[x]
            π_1 = .998340 * π_m    # to account for atmospheric refraction
            eclipse_date = angular_separations[x][1]
            alignment = angular_separations[x][2]

            if alignment == "opposition":
                Total = 1.02 * (π_1 + π_s - s_s) - s_m        # eq 11.13
                Partial = (1.02 * (π_1 + π_s - s_s)) + s_m    # eq 11.14
                Penumbral = (1.02 * (π_1 + π_s + s_s)) + s_m  # eq 11.15
                if sigma < Total:
                    eclipses.append(("Total Lunar Eclipse", eclipse_date, True))
                elif sigma < Partial:
                    eclipses.append(("Partial Lunar Eclipse", eclipse_date,
                                     True))
                elif sigma < Penumbral:
                    eclipses.append(("Penumbral Lunar Eclipse", eclipse_date,
                                     True))
            elif alignment == "conjunction":
                Total = fabs(s_s - s_m + π_m - π_s)           # eq 11.23
                Partial = s_s + s_m + π_m - π_s               # eq 11.21
                if sigma < Total:
                    # if the moon's angular size is smaller than the sun,
                    # the eclipse is annular
                    if s_s > s_m:
                        eclipses.append(("Annular Solar Eclipse",
                                         eclipse_date, False))
                    else:
                        eclipses.append(("Total Solar Eclipse",
                                         eclipse_date, False))
                elif sigma < Partial and less:
                    eclipses.append(("Partial Solar Eclipse",
                                     eclipse_date, False))
    return eclipses


# method angular_semidiameter
# returns the angular semidiameter of an object
# equation 10.47, pg 427, explanatory supplement
# NOTE: equation is wrong in supplement. See erratica
# ----------------
# Inputs: distance (au) of body, novas planet object
# Output: angular semidiameter of the object, in degrees
#
def angular_semidiameter(distance, body):
    if body.number == 10:
        radius = 695700.0   # km
    elif body.number == 11:
        radius = 1738.1     # km
    distance *= const.AU_KM
    return math.degrees(math.asin(radius / distance))


# method parallax
# returns the parallax of an object, which can be very high for the moon
# equation 10.46 / 10.74, page 427 / 434, explanatory supplement
# ----------------
# Input: radius of the object
# Output: horizontal parallax of the object, in degrees
#
def parallax(rp):
    Re = const.ERAD * pow(10, -3)
    rp *= const.AU_KM
    HP = math.asin(Re / fabs(rp))
    return math.degrees(HP)


# method geocentric_ang_sep_finder
# Returns the geocentric angular separation of two objects in opposition
#  or conjunction. See ASTR ALMANAC supplement pg 458, eq 11.9
# ----------------
# Inputs: longitudal motions of moon and sun, relative β_mlination,
#         date, alignment (opposition or conjunction)
# Outputs: A list of geocentric angular separation, date, and alignment
#
def geocentric_ang_sep_finder(SP, SS_prime, βm, date, alignment):
    I_rad = math.radians(I)     # Python math class works in radians
    λ = SP / SS_prime
    I_prime = math.atan((λ / (λ - 1)) * math.tan(I_rad))
    σ = βm * math.cos(I_prime)
    return [σ, date, alignment]




# method m_calculator
# Calculates the angular separation vector between the sun and moon in the
# fundamental shadow plane. x, y, and m are pulled from eq 11.136, pg 492 of
# the explanatory supplement. xdot and ydot I have derived myself.
# ---------------
# Inputs: time, the julian date at a specific instance
# outputs: m and mdot (in degrees), x, and y
def m_calculator(time):
    pos_moon = nov.app_planet(time, moon)
    pos_sun = nov.app_planet(time, sun)
    old_pos_moon = nov.app_planet(time - map_t_step, moon)
    old_pos_sun = nov.app_planet(time - map_t_step, sun)

    # because python's math  works in radians, it is easiest to convert here
    αm = rad((pos_moon[0] * 15.0) + α_correction)
    αs = rad(pos_sun[0] * 15.0)
    δm = rad(pos_moon[1] + δ_correction)
    δs = rad(pos_sun[1])
    a = αs + math.pi
    d = -1 * δs

    # finding the hourly rate of variation at this instant.
    αmdot = (αm - rad((old_pos_moon[0] * 15) + α_correction)) / map_t_step
    αsdot = (αs - rad(old_pos_sun[0] * 15)) / map_t_step
    δmdot = (δm - rad(old_pos_moon[1] + δ_correction))/ map_t_step
    δsdot = (δs - rad(old_pos_sun[1])) / map_t_step
    adot = αsdot
    ddot = -1 * δsdot

    # as per the explanatory supplement
    x = cos(δm) * sin((αm - a))
    y = sin(δm) * cos(d) - cos(δm) * sin(d) * cos(αm - a)

    # splitting terms for line length reasons
    xdot1 = cos(δm) * cos(αm - a) * (αmdot - αsdot)
    xdot2 = -1 * sin(αm - a) * sin(δm) * δmdot
    xdot = xdot1 + xdot2

    ydot1 = δmdot * cos(δm) * cos(d)
    ydot2 = δmdot * sin(d) * cos(αm - a) * sin(δm)
    ydot3 = (αmdot - αsdot) * (sin(d) * sin (αm - a) * cos(δm))
    ydot4 = -1 * ddot * sin(δm) * sin(d)
    ydot5 = -1 * ddot * cos(d) * cos(αm - a) * cos(δm)
    ydot = ydot1 + ydot2 + ydot3 + ydot4 + ydot5

    m = math.sqrt(x ** 2 + y ** 2)      # eq. 11.137
    mdot = (x * xdot + y * ydot) / m    # self derived
    return [math.degrees(m), math.degrees(mdot), x, y]


# method title_generator
# Returns a string of the tital of the eclipse to be plotted
# The choices made here reflect a desire to copy the original fortran plots
# as closely as possible
# Inputs: Event (type of eclipse), lunar (boolean), P1 contact time,
#         P4 contact time
# Outputs: title string, date of the event (year, month, day(s)),
#          output filename
#
def title_generator(event, n, lunar, P1, P4):
    title, filename = '', ''
    n += 1  # to convert from zero-based indexing
    # The following code is messy, couldn't find a Roman numerals library
    if n <= 3:
        for x in range(n):
            title += ("I")
    elif n == 4:
        title += ("IV")
    elif n >= 5:
        title += "V"
        for x in range(5, n):
            title += ("I")
    title += (". - ")
    # The title
    if event.__contains__("Partial"):
        title += ("Partial")
    else:
        if lunar:
            if event.__contains__("Penumbral"):
                title += ("Penumbral")
            else:
                title += ("Total")
        else:
            if event.__contains__("Annular"):
                title += "Annular"
            else:
                title += ("Total")
    title += (" Eclipse of the " )
    if lunar:
        title += ("Moon")
        filename += "L"
    else:
        title += ("Sun")
        filename += "S"
    Greg_date = nov.cal_date(P1)
    year = Greg_date[0]
    month = calendar.month_name[Greg_date[1]]
    day = Greg_date[2]
    plot_date = str(year) + " " + month + " " + str(day)
    filename += str(year) + month[:3] + str(day)

    # This looks to see if the eclipse begins and ends on different days.
    # This is reflected in the title.
    if nov.cal_date(P1)[2] < nov.cal_date(P4)[2]:
        plot_date += "-" + str(Greg_date[2] + 1)

    return title, plot_date, filename


# method RA_time_str
# Returns a string of time of least right ascension
# Inputs: julian time
# Outputs: hour and minute of the conjunction/opposition time,
# to nearest tenth of minute
#
def RA_time_str(jtime):
    G_date = nov.cal_date(jtime)
    month = str(calendar.month_name[G_date[1]]) + " "
    day = str(G_date[2]) + "d "
    hour = math.floor(G_date[3])
    minute = math.floor((G_date[3] - hour) * 60)
    second = round(((G_date[3] - hour) * 60 - minute) * 60, 3)
    pre = "UT of geocentric opposition in RA: "
    string = pre + month + day + str(hour) + "h " + str(minute) +\
           "m " + str(second) + "s"
    return string


# method contact_time_str
# Returns a string of contact times
# Inputs: julian time
# Outputs: hour and minute of the contact time, to nearest tenth of minute
#
def contact_time_str(jtime):
    normal_time = nov.cal_date(jtime)[3]
    hr = math.floor(normal_time)
    minute = round((normal_time - hr) * 60, 1)
    return str(hr) + "h " + str(minute) + "m"


# This class is pulled from a stackexchange answer:
# https://stackoverflow.com/questions/19353576/curved-text-rendering-in-matplotlib
# Courtesy of user Thomas Kühn
# It interacts in much the same way that other pyplot functions do
# for example code, see either mine or follow the above link
#
class CurvedText(mtext.Text):
    """
    A text object that follows an arbitrary curve.
    """
    def __init__(self, x, y, text, axes, **kwargs):
        super(CurvedText, self).__init__(x[0],y[0],' ', axes, **kwargs)

        axes.add_artist(self)

        ##saving the curve:
        self.__x = x
        self.__y = y
        self.__zorder = self.get_zorder()

        ##creating the text objects
        self.__Characters = []
        for c in text:
            if c == ' ':
                ##make this an invisible 'a':
                t = mtext.Text(0,0,'a')
                t.set_alpha(0.0)
            else:
                t = mtext.Text(0,0,c, **kwargs)

            #resetting unnecessary arguments
            t.set_ha('center')
            t.set_rotation(0)
            t.set_zorder(self.__zorder +1)

            self.__Characters.append((c,t))
            axes.add_artist(t)


    ##overloading some member functions, to assure correct functionality
    ##on update
    def set_zorder(self, zorder):
        super(CurvedText, self).set_zorder(zorder)
        self.__zorder = self.get_zorder()
        for c,t in self.__Characters:
            t.set_zorder(self.__zorder+1)

    def draw(self, renderer, *args, **kwargs):
        """
        Overload of the Text.draw() function. Do not do
        do any drawing, but update the positions and rotation
        angles of self.__Characters.
        """
        self.update_positions(renderer)

    def update_positions(self,renderer):
        """
        Update positions and rotations of the individual text elements.
        """

        #preparations

        ##determining the aspect ratio:
        ##from https://stackoverflow.com/a/42014041/2454357

        ##data limits
        xlim = self.axes.get_xlim()
        ylim = self.axes.get_ylim()
        ## Axis size on figure
        figW, figH = self.axes.get_figure().get_size_inches()
        ## Ratio of display units
        _, _, w, h = self.axes.get_position().bounds
        ##final aspect ratio
        aspect = ((figW * w)/(figH * h))*(ylim[1]-ylim[0])/(xlim[1]-xlim[0])

        #points of the curve in figure coordinates:
        x_fig,y_fig = (
            np.array(l) for l in zip(*self.axes.transData.transform([
            (i,j) for i,j in zip(self.__x,self.__y)
            ]))
        )

        #point distances in figure coordinates
        x_fig_dist = (x_fig[1:]-x_fig[:-1])
        y_fig_dist = (y_fig[1:]-y_fig[:-1])
        r_fig_dist = np.sqrt(x_fig_dist**2+y_fig_dist**2)

        #arc length in figure coordinates
        l_fig = np.insert(np.cumsum(r_fig_dist),0,0)

        #angles in figure coordinates
        rads = np.arctan2((y_fig[1:] - y_fig[:-1]),(x_fig[1:] - x_fig[:-1]))
        degs = np.rad2deg(rads)


        rel_pos = 10
        for c,t in self.__Characters:
            #finding the width of c:
            t.set_rotation(0)
            t.set_va('center')
            bbox1  = t.get_window_extent(renderer=renderer)
            w = bbox1.width
            h = bbox1.height

            #ignore all letters that don't fit:
            if rel_pos+w/2 > l_fig[-1]:
                t.set_alpha(0.0)
                rel_pos += w
                continue

            elif c != ' ':
                t.set_alpha(1.0)

            #finding the two data points between which the horizontal
            #center point of the character will be situated
            #left and right indices:
            il = np.where(rel_pos+w/2 >= l_fig)[0][-1]
            ir = np.where(rel_pos+w/2 <= l_fig)[0][0]

            #if we exactly hit a data point:
            if ir == il:
                ir += 1

            #how much of the letter width was needed to find il:
            used = l_fig[il]-rel_pos
            rel_pos = l_fig[il]

            #relative distance between il and ir where the center
            #of the character will be
            fraction = (w/2-used)/r_fig_dist[il]

            ##setting the character position in data coordinates:
            ##interpolate between the two points:
            x = self.__x[il]+fraction*(self.__x[ir]-self.__x[il])
            y = self.__y[il]+fraction*(self.__y[ir]-self.__y[il])

            #getting the offset when setting correct vertical alignment
            #in data coordinates
            t.set_va(self.get_va())
            bbox2  = t.get_window_extent(renderer=renderer)

            bbox1d = self.axes.transData.inverted().transform(bbox1)
            bbox2d = self.axes.transData.inverted().transform(bbox2)
            dr = np.array(bbox2d[0]-bbox1d[0])

            #the rotation/stretch matrix
            rad = rads[il]
            rot_mat = np.array([
                [math.cos(rad), math.sin(rad)*aspect],
                [-math.sin(rad)/aspect, math.cos(rad)]
            ])

            ##computing the offset vector of the rotated character
            drp = np.dot(dr,rot_mat)

            #setting final position and rotation:
            t.set_position(np.array([x,y])+drp)
            t.set_rotation(degs[il])

            t.set_va('center')
            t.set_ha('center')

            #updating rel_pos to right edge of character
            rel_pos += w-used


# method lunar eclipse mapper
# Does all of the calculations to find contact times, positional coordinates,
# etc relevant to the chosen eclipse. Plots it all in matplotlib and
# basemap. This is a long one, so buckle up.
# ----------------
# Inputs: the list of eclipses in the chosen year, the user's
# choice of eclipse number, and the specified ∆T value
#
# Outputs: Three png images -- a map of the moon during the eclipse describing
# its motion, a map of the Earth showing where the eclipse will be visible,
# and a third image combining the previous two -- the complete plot as seen
# online and in the almanac supplement
#
def lunar_eclipse_mapper(eclipses, num, ΔT):
    event = eclipses[num][0]    # string describing the chosen eclipse
    date = eclipses[num][1]     # rough date of eclipse occurence
    start = date - .5           # a start to our search time
    time = start                # the iterative variable
    end = date + .5             # the end to our search time
    least_m_dot_m, least_RA_sep = 100000, 100000    # storage variables
    max_timestep = 12           # the lowest our timestep will go
                                # which is 10^-12 Julian days

    # Classifying the eclipse with easy to reference booleans
    # Setting up contact numbers
    Penumbral, Partial, Total = False, False, False
    if event.__contains__("Partial"):
        patial = True
        n_contacts = [0,1,'MID',5,6]
    elif event.__contains__("Penumbral"):
        Penumbral = True
        n_contacts = [0,'MID', 6]
    else:
        Total = True
        n_contacts = [0,1,2,'MID',4,5,6]

    # iterating to find the midpoint or point of greatest obscuration
    # This loop is derived from section 11.4.2.3, page 493
    # of the explanatory supplement.
    # The smallest m * mdot is saved and this value is calculated at each step
    # as long as it is getting smaller, we continue in that time direction
    # if at any point the calculated m * mdot is larger than
    # the previous calculated one, we switch direction and reduce our timestep
    # the timestep becomes extremely small.
    # Effect on computing time is negligable.
    last_mm = 10000000
    sign = 1
    timestep = .1
    # we cut timestep by a factor of 10 12 times
    for x in range(max_timestep):
        y = 0
        # at most, 100 iterations in one direction for any particular timestep
        while y < 100:
            mdata = m_calculator(time)
            m, mdot = mdata[0], mdata[1]
            mm = fabs(m * mdot)
            if mm * mdot < least_m_dot_m:
                least_m_dot_m = mm
                UT_obsc = time
            # if we're on track, keep going
            if mm < last_mm:
                last_mm = mm
            # if things go wrong, switch and end the y loop.
            elif mm > last_mm:
                sign *= -1
                y = 100
            # increment time regardless
            time = round(time + timestep * sign, x+1)
            y += 1
        timestep /= 10

    ######################
    # The same process for opposition time,
    # based on difference in right ascension
    time = start
    last_diff = 10000000
    sign = 1
    timestep = .1
    for x in range(max_timestep):
        y = 0
        while y < 100:
            alpha_s = nov.app_planet(time, sun)[0] * 15
            alpha_m = (nov.app_planet(time, moon)[0] * 15) + α_correction
            RA_diff = fabs(alpha_m - alpha_s) - 180

            if RA_diff < least_RA_sep:
                least_RA_sep = RA_diff
                UT_opp = time
            if RA_diff < last_diff:
                last_diff = RA_diff
            elif RA_diff > last_diff:
                sign *= -1
                y = 100
            time = round(time + timestep * sign, x+1)
            y += 1
        timestep /= 10

    # Saving the values and applying ∆T
    UT_obsc -= ΔT
    UT_opp -= ΔT

    # Finally, we find the contact times of specific events
    # The math described here is found in section 11.4.2.2, page 493 of the
    # explanatory supplement. It is relatively simple but requires a specific
    # choice of time guesses to be mathematically viable.
    contact_times = {}
    for n in n_contacts:
        # if the contact occurs before midpoint, we start our time iteration
        # at the far front end of the eclipse - well before it happens
        if type(n) == int and n < 4:
            time = start
        elif n == 'MID':
            time= UT_obsc
        # Otherwise, we start it after it is done
        else:
            time = end

        # assigning strings to the contact points
        if n == 0 or n == 6:
            c = "P1"
            if c in contact_times:
                c = "P4"
        elif n == 1 or n == 5:
            c = "U1"
            if c in contact_times:
                c = "U4"
        elif n == 2 or n == 4:
            c = "U2"
            if c in contact_times:
                c = "U3"

        # At each step, the size and parallax of each body are determined
        # for maximum accuracy -- because we have a computer, so why not?
        # One of the two roots of the quadratic described -- the mathematically
        # viable one which has a chance of approaching 0 -- is taken
        # Because L and m0 are both absolute values, they can never be negative
        # thus the other root t = (L + m0) / m0dot never provides a real answer
        for x in range(10):
            # m and mdot are calculated
            mdata = m_calculator(time)
            m0, m0dot = mdata[0], mdata[1]
            # the umbral and penumbral shadow sizes are found
            moon_d = nov.app_planet(time, moon)[2]
            sun_d = nov.app_planet(time, sun)[2]
            s_m = angular_semidiameter(moon_d, moon)
            s_s = angular_semidiameter(sun_d, sun)
            pi_m = parallax(moon_d)
            pi_s = parallax(sun_d)
            pi_1 = .998340 * pi_m
            f1 = 1.02 * (pi_1 + pi_s + s_s)
            f2 = 1.02 * (pi_1 + pi_s - s_s)
            # value of L is calculated based on contact point
            if n == 0 or n == 6:
                L = f1 + s_m
            elif n == 1 or n == 5:
                L = f2 + s_m
            elif n == 2 or n == 4:
                L = f2 - s_m
            # we take a break from the usual loop and sneak in a quick calc
            # of the magnitude of the eclipse for our plot
            # formula here is found in section 11.4.2.4, pg 494
            # of the explanatory supplement
            elif n == 'MID':
                if Penumbral:
                    Li = f1 + s_m
                else:
                    Li = f2 + s_m
                mag = round((Li - m0) / (2 * s_m), 3)
                contact_times[n] = time
                break   # no need to carry on any calculations for the midpoint

            # t is calculated, from this we get a new time, the loop repeats
            # this typically converges within 4-5 iterations.
            # We do 10 for safety reasons
            t = (L - m0) / m0dot
            time += t

        # finally, we save the contact point and time to a dictionary
        # after applying ∆T
        if n is not 'MID':
            contact_times[c] = time - ΔT



    ####################
    #     PLOTTING     #
    ####################

    print(" Plotting...")
    print(" Please ignore deprecation warnings")
    # Finding the sizes of the umbral and penumbral shadows
    # at the time of greatest obscuration
    # for obvious reasons, they can only have one value from here on out
    moon_d = nov.app_planet(UT_obsc, moon)[2]
    sun_d = nov.app_planet(UT_obsc, sun)[2]
    s_m_plotting = angular_semidiameter(moon_d, moon)
    s_s = angular_semidiameter(sun_d, sun)
    pi_m = parallax(moon_d)
    pi_s = parallax(sun_d)
    pi_1 = .998340 * pi_m
    f1 = 1.02 * (pi_1 + pi_s + s_s)

    # drawing the penumbral and umbral circles
    c1 = plt.Circle((0, 0), radius=f1, color=".9", zorder=1)
    c2 = plt.Circle((0, 0), radius=f2, color=".65", zorder=1)

    fig, ax = plt.subplots()
    ax.add_artist(c1)
    ax.add_artist(c2)

    # axes data
    plt.axis("scaled")
    plt.axis("off")
    # setting the axes limits
    # has to be done after scaling
    plt.xlim(-f1 - 1, f1 + 1)
    plt.ylim(-f1 - 1, f1 + 1)

    # Labeling penumbral and umbral circles
    # zorder must be set manually
    # spacing in string is a necessary evil
    # curves that are traced are declared in-function
    text = CurvedText(x = -f2 * np.cos(np.linspace(0,2*np.pi,100)),
                      y = f2 * np.sin(np.linspace(0,2*np.pi,100)),
                      text="           UMBRA", fontname = "Times New Roman",
                      va = "top", ha='center', axes=ax, fontsize=6)
    text.set_zorder(2)
    text = CurvedText(x = -f1 * np.cos(np.linspace(0,2*np.pi,100)),
                      y = f1 * np.sin(np.linspace(0,2*np.pi,100)),
                      text="                    PENUMBRA",
                      fontname = "Times New Roman",
                      va = "top", axes=ax, fontsize=6)
    text.set_zorder(2)
    text = CurvedText(x = -f2 * np.cos(np.linspace(0,2*np.pi,100)),
                      y = -f2 * np.sin(np.linspace(0,2*np.pi,100)),
                      text="          UMBRA", fontname = "Times New Roman",
                      va = "bottom", ha='center', axes=ax, fontsize=6)
    text.set_zorder(2)
    text = CurvedText(x = -f1 * np.cos(np.linspace(0,2*np.pi,100)),
                      y = -f1 * np.sin(np.linspace(0,2*np.pi,100)),
                      text="                    PENUMBRA",
                      fontname = "Times New Roman",
                      va = "bottom", ha='center', axes=ax, fontsize=6)
    text.set_zorder(2)

    # plotting the central cross
    ax.annotate(xy=(0, 0), xycoords='data', s="+", color='blue',
                ha='center', va='center')



    # Drawing the title and date
    RA_time_label = RA_time_str(UT_opp)
    title, date, filename = title_generator(event, num, True,
                              contact_times["P1"], contact_times["P4"])
    ax.annotate(xy=(0, 1), textcoords='axes fraction', s=title, fontsize=9.5,
                fontweight='bold', fontname='Times New Roman',
                ha='left', va='bottom')
    ax.annotate(xy=(1, 1), textcoords='axes fraction', s=date, fontsize=9.5,
                fontweight='bold', fontname='Times New Roman',
                ha='right', va='bottom')


    # Some positional information for our graphics
    β_m_sum = 0             # to find avg latitude of moon
    moon_on_top = False
    for contact in contact_times:
        delta_m = nov.app_planet(contact_times[contact], moon)[1]
        delta_s = nov.app_planet(contact_times[contact], sun)[1]
        β_m = delta_s + delta_m
        if contact == 'P4':
            # placing the compass away from the moon on the left
            if β_m < 0:
                compass_xy = (.05, .7)
            else:
                compass_xy = (.05, .3)
        if contact == 'P1':
            # The same idea, for the 30 arcsec legend
            if β_m < 0:
                legend_xy = (f1+.828, f1+.25)
                linex, liney = [f1+.8, f1+.8], [f1, f1+.5]
            else:
                legend_xy = (f1+.828, -f1)
                linex, liney = [f1+.8, f1+.8], [-f1-.25, -f1+.25]
        β_m_sum += β_m
    if β_m_sum / contact_times.values().__len__() > 0:
        moon_on_top = True

    # Plotting the compass
    fn = get_sample_data(compass_img_file, asfileobj=False)
    arr_img = plt.imread(fn, format='png')
    imagebox = OffsetImage(arr_img, zoom=0.1)
    imagebox.image.axes = ax
    ab = AnnotationBbox(imagebox, xy=compass_xy, xycoords='axes fraction',
                        frameon=False)
    ax.add_artist(ab)

    # Plotting the legend line
    ax.plot(linex, liney, linewidth=1, color='k')
    ax.annotate(xy=legend_xy, s="30 arc-minutes", rotation=270, fontsize=4,
                ha='center', va='center')

    # Annotating the UT opposition and magnitude strings, if moon is on top
    # If not, they will be plotted with the world map
    if Penumbral:
        mag_str = "Penumbral"
    else:
        mag_str = "Umbral"
    mag_str += " magnitude of the eclipse: " + str(mag)
    if moon_on_top:
        ax.annotate(xy=(0, .98), textcoords='axes fraction', s=RA_time_label,
                    fontsize=7, fontname='Times New Roman',
                    ha='left', va='bottom')
        ax.annotate(xy=(1, .98), textcoords='axes fraction', s=mag_str,
                    fontsize=7, fontname='Times New Roman',
                    ha='right', va='bottom')

        sign = 1
        for x in range(contact_times.__len__() - 1):
                offset = 0
                contact = list(contact_times.keys())[x]
                t = list(contact_times.values())[x]
                new_t = list(contact_times.values())[x+1]
                ra = (nov.app_planet(t, moon)[0] * 15) + α_correction
                ra1 = (nov.app_planet(new_t, moon)[0] * 15) + α_correction

                if fabs(ra1 - ra) < .3:
                    contact_times[contact] = [t, -.2]
                    contact_times['MID'] = [new_t, 0]
                    contact_times[list(contact_times.keys())[x+2]] = [list(contact_times.values())[x+2], .2]
                    break


    ################################################
    # Plotting the actual moon at each contact point
    opp_alpha = (nov.app_planet(UT_opp, moon)[0] * 15) + α_correction
    alphas, deltas = [], []
    # Loading moon image
    fn = get_sample_data(moon_img_file, asfileobj=False)
    fn_red = get_sample_data(red_moon_img_file, asfileobj=False)
    arr_moon = plt.imread(fn, format='png')
    arr_moon_red = plt.imread(fn_red, format='png')
    imsize = (s_m_plotting / f1) * .12

    for contact in contact_times:
        t = contact_times[contact]
        offset = 0

        # drawing ephemeris data for the contact points
        try:
            moon_data = nov.app_planet(t, moon)
            sun_data = nov.app_planet(t, sun)
        except TypeError:
            t = contact_times[contact][0]
            offset = contact_times[contact][1]
            moon_data = nov.app_planet(t, moon)
            sun_data = nov.app_planet(t, sun)

        lineoffset = 0
        if contact == 'MID':
            lineoffset = .2

        # Specifying data
        alpha_m = (moon_data[0] * 15) + α_correction
        delta_m = moon_data[1] + δ_correction
        dist_m = moon_data[2]
        alpha_s, delta_s = sun_data[0] * 15, sun_data[1]

        # finding coordinates of the moon at contact point
        β_m = delta_s + delta_m
        RA_diff = alpha_m - alpha_s
        if RA_diff < 0:
             RA_diff += 180
        else:
            RA_diff -= 180
        RA_diff *= -1
        # saving them for plotting the red path line later
        alphas.append(RA_diff)
        deltas.append(β_m)

        # Position Angle calculation
        mdata = m_calculator(t)
        x, y = mdata[2], mdata[3]
        PA = round(math.degrees(math.atan(x/y)), 1)
        if PA < 0:
            PA += 360
        elif PA > 360:
            PA -= 360

        # Building the constant string for annotation
        contact_str = contact + "\n" + contact_time_str(t) + "\n"
        mark_type = ""      # the red crosses and circles
        cmap = arr_moon
        if contact is "P1" or contact is "P4":
            if Penumbral:
                contact_str += "P.A. " + str(PA)
            mark_type = "+"
        elif contact is "U1" or contact is "U4":
            contact_str += "P.A. " + str(PA + 180)
            mark_type = "☉"
        elif contact is "U2" or contact is "U3":
            mark_type = "⊕"
            if Total:
                cmap = arr_moon_red
        elif contact is "MID" and Total:
            cmap = arr_moon_red

        # Saving the mark for later, when we plot onto the world map
        contact_times[contact] = [t, mark_type, offset]
        # contact string position is dependent on moon position
        if moon_on_top:
            contact_line_coordsx = [RA_diff, RA_diff]
            contact_line_coordsy = [β_m, -f1 - .5 + fabs(offset) + lineoffset]
            contact_mark_coords = (RA_diff, -f1 - .25)
            contact_str_coords = (RA_diff, -f1 - .55 + .025 + offset)
            contact_vas = 'top'
        else:
            contact_line_coordsx= [RA_diff, RA_diff]
            contact_line_coordsy = [β_m, f1 + .25 - .01 - fabs(offset) - lineoffset]
            contact_mark_coords = (RA_diff, f1 + .15)
            contact_str_coords = (RA_diff, f1 + .25 + offset)
            contact_vas = 'bottom'

        # Plotting the dotted line
        ax.plot(contact_line_coordsx, contact_line_coordsy, linestyle=":",
                color="red", linewidth=.5, zorder=2)
        # Annotating the mark / contact symbol
        ax.annotate(mark_type, xy=contact_mark_coords, ha='center', va='center',
                    fontsize=10, color='red', zorder=2)
        # Annotating the actual contact string
        ax.annotate(contact_str, xy=contact_str_coords, ha='center',
                    va=contact_vas, fontsize=6, fontname="Times New Roman",
                    zorder=3)


        # Finally, plotting the image of the moon
        imbox = OffsetImage(arr=cmap, zoom=imsize)
        imbox.image.axes = ax
        ab = AnnotationBbox(imbox, xy=(RA_diff, β_m), xycoords='data',
                            frameon=False)
        ax.add_artist(ab)

        # end contact loop

    # Calculating the slope and y-intercept of the path of the moon
    rise = deltas[-1] - deltas[0]
    run = alphas[-1] - alphas[0]
    slope = rise / run
    b = deltas[0] - slope * alphas[0]

    # plotting the red arrows denoting moon motion
    arrow1x = alphas[0] + .5
    arrow1y = slope * arrow1x + b
    arrow2x = alphas[-1] - .5
    arrow2y = slope * arrow2x + b
    ax.annotate(xy=(arrow1x, arrow1y), s="►",
                rotation=math.degrees(math.atan(slope)) - 180,
                fontsize=10, ha='center', va='center', color='red')
    ax.annotate(xy=(arrow2x, arrow2y), s="►",
                rotation=math.degrees(math.atan(slope)) - 180,
                fontsize=10, ha='center', va='center', color='red')

    # plotting the red moon path
    for x in range(-5, 5):
        deltas.append(slope * x + b)
        alphas.append(x)
    ax.plot(alphas, deltas, "-r", linewidth=.5)

    # plotting the ecliptic
    slope = math.degrees(math.atan(slope))
    start_δ = fabs(nov.app_planet(contact_times['P1'][0], moon)[1])
    end_δ = fabs(nov.app_planet(contact_times['P4'][0], moon)[1])
    if end_δ < start_δ and end_δ > 0 and start_δ > 0:
        ecl_slope = math.tan(rad(slope - I))
    else:
        ecl_slope = math.tan(rad(slope + I))

    ecl_x, ecl_y = [], []
    for x in range(-5, 5):
        ecl_x.append(x)
        ecl_y.append(ecl_slope * x)
    ax.plot(ecl_x, ecl_y, '-b', linewidth=.5, zorder=3)
    e_x1, e_x2 = -f1-.5, f1+.5
    e_y1, e_y2 = e_x1 * ecl_slope, e_x2 * ecl_slope
    e_y_offset = -.04
    if moon_on_top:
        e_y_offset = .04
    # ecliptic labels
    ax.annotate(xy=(e_x1, e_y1 + e_y_offset), s="Ecliptic",
                rotation=math.degrees(math.atan(ecl_slope)),
                ha='center', va='center', color='blue',
                fontname="Times New Roman", fontsize=3.5, zorder=3)
    ax.annotate(xy=(e_x2, e_y2 + e_y_offset), s="Ecliptic",
                rotation=math.degrees(math.atan(ecl_slope)),
                ha='center', va='center', color='blue',
                fontname="Times New Roman", fontsize=3.5, zorder=3)

    # resizing image, saving the figure, clearing out pyplot
    plt.subplots_adjust(left=.05, bottom=0, right=.95, top=.95,
                        wspace=0, hspace=0)
    plt.savefig(filename + "_moonmap.png", dpi=650)
    plt.clf()


    #############################
    #    Plotting the Earth     #
    #############################

    # finding central coordinates for centering the map
    midt = contact_times['MID'][0]
    center_data = nov.app_planet(midt, moon)
    φ_mid = center_data[1]
    λ_mid = (nov.sidereal_time(int(midt), midt%1, ΔT) - center_data[0]) * -15
    if λ_mid < 0:
        λ_mid += 360
    center_data = nov.app_planet(midt, moon)

    # Drawing the map, borders, meridians and parallels
    earth = Basemap(projection='cyl', lon_0=λ_mid, lat_0=φ_mid,
                    llcrnrlon=λ_mid-180, llcrnrlat=-80,
                    urcrnrlat=80, urcrnrlon =λ_mid+180,
                    suppress_ticks=True)
    earth.drawcoastlines(color='blue', linewidth=.25)
    earth.drawcountries(color='green', linewidth=.25)
    earth.drawmapboundary()
    parallels = np.arange(-100, 100 , 20)
    meridians = np.arange(-180, 180, 30)
    earth.drawparallels(parallels, dashes=[3, 409],
                        labels=[1, 1, 1, 1], fontsize=7)
    earth.drawmeridians(meridians, dashes=[3, 178],
                        labels=[1, 1, 1, 1], fontsize=7)


    longs, lats = [], []       # for red dotted line between contact points
    dark_top = True            # for determing label placement after loop
    del contact_times["MID"]   # Mid is no longer necessary
    maxcolor = .5 + (.075 * contact_times.__len__())
    labelcount = 0             # for iterating horizon labels up or down
    for contact in contact_times:
        t = contact_times[contact][0]
        mark_type = contact_times[contact][1]
        offset = contact_times[contact][2]
        index = list(contact_times.values()).index([t, mark_type, offset])

        # Finding the two colors of gray we need
        color1 = .5 + (.075 * index)
        color2 = maxcolor - (.075 * index) - .075
        color1, color2 = str(color1), str(color2)

        # if the moon is high, shading appears on the bottom of the plot
        # and the horizon labels appear from north to south vs south to north
        sign = 1
        if nov.app_planet(t, moon)[1] > 0:
            sign = -1
            dark_top = False

        # finding the coordinates at which the moon is at zenith, as described
        # in section 11.4.2.6 pg 494 of the explanatory supplement
        # NOTE: longitude is west-positive (ugh)
        # NOTE: sidetime requires those inputs to correctly return Greenwich
        # apparent sidereal time -- do not touch it!!!
        # NOTE: as sidetime is 0-24 hours, coordinates can sometimes appear
        # in seperate 360˚ domains. This is corrected manually below
        sidetime = nov.sidereal_time(int(t), t%1, ΔT)
        αm = nov.app_planet(t, moon)[0]
        λ0 = 15 * (sidetime - αm) - α_correction
        φ0 = nov.app_planet(t, moon)[1]
        λ0 *= -1
        if λ0 > (λ_mid + 180):
            λ0 -= 360
        elif λ0 < (λ_mid - 180):
            λ0 += 360
        plt.annotate(xy=(λ0, φ0), s=mark_type, fontsize=6, color='red',
                     ha='center', va='center')
        # we save these to plot the red dotted line
        longs.append(λ0), lats.append(φ0)

        # The horizons as described by section 11.4.2.6, pg 494 of the
        # explanatory supplement. Each 'half' is plotted seperately, west
        # and east. This simplifies the following code immensly.
        circle_longs1 = np.arange(λ0, λ0 + 250)
        circle_longs2 = np.arange(λ0 - 250, λ0)
        great_circle_φ1, great_circle_φ2 = [], []

        # Calculating the label coordinates, using both number of contacts
        # and positions of moon zenith points to determine where on the map
        # they will appear. If the moon appears more to the northern hemisphere,
        # the labels will start at 20S and rise to 20N
        # labels are offset by 5˚ longitude for aesthetic effect
        label_lat1 = sign * (20 - 40 / (contact_times.__len__() - 1) * labelcount)
        label_long1 = λ0 - math.degrees(math.acos(-1 * math.tan(rad(label_lat1)) * math.tan(rad(φ0))))
        label_lat2 = sign * (-20 + 40 / (contact_times.__len__() - 1) * labelcount)
        label_long2 = (λ0) + math.degrees(math.acos(-1 * math.tan(rad(label_lat2)) * math.tan(rad(φ0))))
        label_xy1 = (label_long1 + 5, label_lat1)
        label_xy2 = (label_long2 - 5, label_lat2)
        labelcount += 1

        # plotting the great circle -- east, then doing the shading with the
        # appropriate color in the appropriate direction
        for step in circle_longs1:
            φ = math.degrees(math.atan(-1 * (1.0 / math.tan(rad(φ0)) * cos(rad(λ0 - step)))))
            great_circle_φ1.append(φ)
        plt.plot(circle_longs1, great_circle_φ1, linewidth=.5, color='red')
        plt.fill_between(circle_longs1, great_circle_φ1, 90 * sign, color=color1, zorder = 0 - index)

        # same but to the west
        for step in circle_longs2:
            φ = math.degrees(math.atan(-1 * (1.0 / math.tan(rad(φ0)) * cos(rad(λ0 - step)))))
            great_circle_φ2.append(φ)
        plt.plot(circle_longs2, great_circle_φ2, linewidth=.5, color='red')
        plt.fill_between(circle_longs2, great_circle_φ2, 90 * sign, color=color2)

        # finally comes the label placement
        plt.annotate(s=contact, xy=label_xy1, xycoords='data',
                     va='center', ha='center', fontname='Times New Roman',
                     alpha=1, fontsize=6,
                     bbox=dict(boxstyle="square", fc="w", pad=0, ec='w'))
        plt.annotate(s=contact, xy=label_xy2, xycoords='data',
                     va='center', ha='center', fontname='Times New Roman',
                     alpha=1, fontsize=6,
                     bbox=dict(boxstyle="square", fc="w", pad=0, ec='w'))

    # plotting the red dotted line
    plt.plot(longs, lats, linestyle="dashed", color='red', linewidth=.25)

    # Plotting labels denoting moon visibility
    if dark_top:
        no_eclipse_vis_xy1 = (.1, .95)
        no_eclipse_vis_xy2 = (.9, .95)
        risingxy = (.25, .05)
        settingxy = (.75, .05)
        fullxy = (.5, .05)
    else:
        no_eclipse_vis_xy1 = (.1, .1)
        no_eclipse_vis_xy2 = (.9, .1)
        risingxy = (.25, .9)
        settingxy = (.75, .9)
        fullxy = (.5, .9)
    plt.annotate(s='Full eclipse visible', xy=fullxy, xycoords='axes fraction',
                 va='center', ha='center', fontname='Times New Roman',
                 alpha=1, fontsize=6,
                 bbox=dict(boxstyle="square", fc="w", pad=.2, ec='w'))
    plt.annotate(s='No eclipse visible', xy=no_eclipse_vis_xy1,
                 xycoords='axes fraction', va='center', ha='center',
                 fontname='Times New Roman', alpha=1, fontsize=6,
                 bbox=dict(boxstyle="square", fc="w", pad=.2, ec='w'))
    plt.annotate(s='No eclipse visible', xy=no_eclipse_vis_xy2,
                 xycoords='axes fraction', va='center', ha='center',
                 fontname='Times New Roman', alpha=1, fontsize=6,
                 bbox=dict(boxstyle="square", fc="w", pad=.2, ec='w'))
    plt.annotate(s='Moon rising during eclipse', xy=risingxy,
                 xycoords='axes fraction', va='center', ha='center',
                 fontname='Times New Roman', alpha=1, fontsize=6,
                 bbox=dict(boxstyle="square", fc="w", pad=.2, ec='w'))
    plt.annotate(s='Moon setting during eclipse', xy=settingxy,
                 xycoords='axes fraction', va='center', ha='center',
                 fontname='Times New Roman', alpha=1, fontsize=6,
                 bbox=dict(boxstyle="square", fc="w", pad=.2, ec='w'))

    # Creating the ∆T string for labeling
    dt_str = "∆T = "
    if ΔT > 0:
        dt_str += "+"
    else:
        dt_str += "-"
    dt_str += str(round(ΔT * 60 * 60 * 24, 1)) + "s"
    # Annotating miscellanious labels
    # ∆T, ©HM Nautical Almanac Office, worldmap title
    plt.annotate(xy=(0, -.1), xycoords='axes fraction',
                 s="©HM Nautical Almanac Office",
                 fontname="Times New Roman", fontsize=6,
                 ha='center', va='bottom')
    plt.annotate(xy=(.5, -.1), xycoords='axes fraction',
                 s="Areas of visibility of the eclipse at different stages",
                 fontsize=8, fontname="Times New Roman",
                 ha='center', va='bottom')
    plt.annotate(xy=(1, -.1), xycoords='axes fraction', s=dt_str,
                 fontsize=6, fontname="Times New Roman",
                 ha='center', va='bottom')

    # Plotting the RA string and magnitude, if it was not done earlier,
    # dependent on moon position
    if moon_on_top == False:
        plt.annotate(xy=(0, 1.1),xycoords='axes fraction',
                     s=RA_time_label, fontsize=7, fontname='Times New Roman',
                     ha='left', va='center')
        plt.annotate(xy=(1, 1.1), xycoords='axes fraction', s=mag_str,
                     fontsize=7, fontname='Times New Roman',
                     ha='right', va='center')

    # Saving the worldmap
    plt.savefig(filename + "_worldmap.png", dpi=500, bbox_inches='tight')


    # Now that we have a moon map and world map, we marry the two together
    # and do some minor aesthetic editing.
    # Specifically, we fill in blackspace with white pixels, and draw a black
    # border around the entire new image
    im1 = Image.open(filename + "_moonmap.png")
    im2 = Image.open(filename + "_worldmap.png")
    new = Image.new('RGB', (im2.width, im1.height + im2.height))
    new.paste(im1, (int((im2.width-im1.width)/2), 0))
    new.paste(im2, (0, int(im1.height)))
    pixelsNew = new.load()

    for pixelx in range(int((im2.width-im1.width)/2)):
        for pixely in range(im1.height):
            pixelsNew[pixelx, pixely] = (256, 256, 256)

    for pixelx in range(int(im2.width - (im2.width-im1.width)/2), im2.width):
        for pixely in range(im1.height):
            pixelsNew[pixelx, pixely] = (256, 256, 256)

    for pixelx in range(10, 20):
        for pixely in range(20, new.height - 20):
            pixelsNew[pixelx, pixely] = (0, 0, 0)

    for pixelx in range(new.width - 20, new.width - 10):
        for pixely in range(20, new.height - 20):
            pixelsNew[pixelx, pixely] = (0, 0, 0)

    for pixelx in range(20, new.width - 20):
        for pixely in range(20, 30):
            pixelsNew[pixelx, pixely] = (0, 0, 0)

    for pixelx in range(20, new.width - 20):
        for pixely in range(new.height - 30, new.height - 20):
            pixelsNew[pixelx, pixely] = (0, 0, 0)

    new.save(filename + ".png")


    #TODO
    # 1. Ecliptic is sometimes + I, sometimes - I, and I'm not sure how to descriminate between the two
    #    not not not fixed damn this is hard
    # 6. Make the umbral moons red
    # size moons based on umbral shadows


# Execute main
if __name__ == '__main__':
    main()
    print("\n\n\n Program terminates normally")
