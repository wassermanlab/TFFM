"""
    :platform: Unix
    :synopsis: Define a lot of functions necessary to draw the logos.

"""


from constants import BLACK, WHITE, BLUE, RED, GREEN, YELLOW, ALPHABET
from constants import TFFM_KIND, LOGO_TYPE
import tffm_module


def print_header(output, height, width):
    """
    Print the header of the svg to the output.

    :arg output: Stream where to write.
    :type output: file
    :arg height: Height of the final svg.
    :type height: int
    :arg width: Width of the final svg.
    :type width: int

    """

    dtd_url = "http://www.w3.org/TR/2001/REC-SVG-20010904/DTD/svg10.dtd"
    svg_url = "http://www.w3.org/2000/svg"
    xlink_url = "http://www.w3.org/1999/xlink"
    output.write("""<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
                 <!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 1.0//EN" "%s">
                 <svg height="%d" width="%d" xmlns="%s" xmlns:xlink="%s">
                 """ % (dtd_url, height, width, svg_url, xlink_url))


def print_footer(output):
    """
    Print the footer of the svg to the output.

    :arg output: Stream where to write.
    :type output: file

    """

    output.write("\n</svg>\n")


def print_ic(output, length, information_content, x_coord, y_coord):
    """
    Print the information content underneath the logo.

    :arg output: Stream where to write.
    :type output: file
    :arg length: Number of positions in the logo.
    :type length: int
    :arg information_content: Value of the information content.
    :type information_content: float
    :arg x_coord: x-coordinate where to write the information content value in
        the svg.
    :type x_coord: float
    :arg y_coord: y-coordinate where to write the information content value in
        the svg.
    :type y_coord: float

    """

    output.write("<text x=\"%f\" y=\"%d\">\n" % (
        x_coord / 2 + 30 * (length - 1), y_coord + 45))
    output.write("IC = %.2f\n</text>\n" % (information_content))


def draw_axis(output, length, x_coord, y_coord):
    """
    Print the svg instructions to draw the axis to the output.

    :arg output: Stream where to write.
    :type output: file
    :arg length: Number of positions in the logo.
    :type length: float
    :arg x_coord: x-coordinate where to write the axis in the svg.
    :type x_coord: float
    :arg y_coord: y-coordinate where to write the axis in the svg.
    :type y_coord: float

    """

    output.write("<line x1=\"%d\" y1=\"%d\" " % (x_coord + 5, y_coord))
    output.write("x2=\"%d\" y2=\"%d\" " % (
        x_coord + 7 + 60 * length, y_coord))
    output.write("style=\"stroke:rgb(99,99,99);stroke-width:2\" />\n")
    for i in xrange(length + 1):
        draw_rectangle(output, x_coord + 5 + 60 * i, y_coord, 7, 2)
        if i < length:
            output.write("<text x=\"%f\" y=\"%d\">\n" % (x_coord + 27 + 60 * i,
                                                         y_coord + 20))
            output.write("%d\n</text>\n" % (i + 1))


def draw_y_axis(output, x_coord=80, y_coord=120):
    """
    Print the svg instructions to draw the axis to the output.

    :arg output: Stream where to write.
    :type output: file
    :arg x_coord: x-coordinate where to start (bottom) the Y-axis in the svg.
    :type x_coord: float
    :arg y_coord: y-coordinate where to start (bottom) the Y-axis in the svg.
    :type y_coord: float

    """

    output.write("<line x1=\"%d\" y1=\"%d\" " % (x_coord, y_coord - 8))
    output.write("x2=\"%d\" y2=\"%d\" " % (x_coord, y_coord - 110))
    output.write("style=\"stroke:rgb(99,99,99);stroke-width:2\" />\n")
    draw_rectangle(output, x_coord - 7, y_coord - 10, 2, 7)
    output.write("<text x=\"%d\" y=\"%d\">\n0\n</text>\n" % (x_coord - 17,
                                                             y_coord - 5))
    draw_rectangle(output, x_coord - 7, y_coord - 60, 2, 7)
    output.write("<text x=\"%d\" y=\"%d\">\n1\n</text>\n" % (x_coord - 17,
                                                             y_coord - 55))
    draw_rectangle(output, x_coord - 7, y_coord - 110, 2, 7)
    output.write("<text x=\"%d\" y=\"%d\">\n2\n</text>\n" % (x_coord - 17,
                                                             y_coord - 105))
    output.write("<text x=\"%d\" y=\"%d\"\n" % (x_coord - 30, y_coord - 55))
    output.write("transform=\"rotate(-90, %d, %d)\">\n" % (x_coord - 30,
                                                           y_coord - 50))
    output.write("bits\n</text>\n")


def draw_ellipse(output, x_center, y_center, x_radius, y_radius, colour=BLACK,
                 opacity=1.0):
    """
    Print the svg instructions to draw an ellipse to the output.

    :arg output: Stream where to write.
    :type output: file
    :arg x_center: x-coordinate of the center of the ellipse.
    :type x_center: float
    :arg y_center: y-coordinate of the center of the ellipse.
    :type y_center: float
    :arg x_radius: Radius of the ellipse on the x-axis.
    :type x_radius: float
    :arg y_radius: Radius of the ellipse on the y-ayis.
    :type y_radius: float
    :arg colour: Color of the ellipse (default: black). The colour is given in
        RGB.
    :type colour: tuple (int, int, int)
    :arg opacity: Opacity of the ellipse (default: 1.0)
    :type opacity: float

    :note: Look at the ellipse statement of an svg file for more information.

    """

    output.write("<ellipse cx=\"%f\" cy=\"%f\" rx=\"%f\" ry=\"%f\" " %
                 (x_center, y_center, x_radius, y_radius))
    output.write("fill=\"rgb%s\" style=\"fill-opacity:%.2f\" />\n" %
                 (str(colour), opacity))


def draw_polygon(output, pts, colour=BLACK, opacity=1.0):
    """
    Print the svg instructions to draw an polygon to the output.

    :arg output: Stream where to write.
    :type output: file
    :arg pts: Coordinates of the points defining the polygon. For each point in
        the list, its coordinates are given by the couple
        (x-coordinate, y-coordinate).
    :type pts: list of tuple (float, float)
    :arg colour: Color of the polygon (default: black). The colour is given in
        RGB.
    :type colour: tuple (int, int, int)
    :arg opacity: Opacity of the polygon (default: 1.0)
    :type opacity: float

    :note: Look at the polygon statement of an svg file for more information.

    """

    output.write("<polygon points=\" ")
    for i in xrange(len(pts) - 1):
        output.write("%f, %f, " % (pts[i][0], pts[i][1]))
    output.write("%f, %f\" " % (pts[len(pts) - 1][0], pts[len(pts) - 1][1]))
    output.write("fill=\"rgb%s\" style=\"fill-opacity:%.2f\" />\n" %
                 (str(colour), opacity))


def draw_rectangle(output, x_coord, y_coord, height, width, colour=BLACK,
                   opacity=1.0):
    """
    Print the svg instructions to draw an rectangle to the output.

    :arg output: Stream where to write.
    :type output: file
    :arg x_coord: x-coordinate of the bottom-left corner of the rectangle.
    :type x_coord: float
    :arg y_coord: y-coordinate of the bottom-left corner of the rectangle.
    :type y_coord: float
    :arg height: Height of the rectangle.
    :type height: float
    :arg width: Width of the rectangle.
    :type width: float
    :arg colour: Color of the rectangle (default: black). The colour is given
        in RGB.
    :type colour: tuple (int, int, int)
    :arg opacity: Opacity of the rectangle (default: 1.0)
    :type opacity: float

    :note: Look at the rectangle statement of an svg file for more information.

    """

    output.write("<rect height=\"%f\" width=\"%f\" " % (height, width))
    output.write("x=\"%f\" y=\"%f\" " % (x_coord, y_coord))
    output.write("fill=\"rgb%s\" style=\"fill-opacity:%.2f\" />\n" %
                 (str(colour), opacity))


def draw_letter_a(output, x_coord, y_coord, width, height, colour=RED,
                  opacity=1.0):
    """
    Print the svg instructions to draw the letter 'A' to the output.

    :arg output: Stream where to write.
    :type output: file
    :arg x_coord: x-coordinate of the bottom left corner of the box containing
        the 'A'.
    :type x_coord: float
    :arg y_coord: y-coordinate of the bottom left corner ot the box containing
        the 'A'.
    :type y_coord: float
    :arg width: Width of the box containing the 'A'.
    :type width: float
    :arg height: Height of the box containing the 'A'.
    :type height: float
    :arg colour: Color of the 'A' (default: black). The colour is given in
        RGB.
    :type colour: tuple (int, int, int)
    :arg opacity: Opacity of the 'A' (default: 1.0)
    :type opacity: float

    :note: The 'A' is drawn by creating a red polygon giving the shape and a
        white internal triangle above to make the bar in the middle of the 'A'
        appear.

    """

    output.write("<!-- Begin 'A' -->\n")
    pt1 = (x_coord, y_coord + height)
    pt2 = (x_coord + width * 0.42, y_coord)
    pt3 = (x_coord + width * 0.58, y_coord)
    pt4 = (x_coord + width, y_coord + height)
    pt5 = (x_coord + 0.85 * width, y_coord + height)
    pt6 = (x_coord + 0.725 * width, y_coord + 0.75 * height)
    pt7 = (x_coord + 0.275 * width, y_coord + 0.75 * height)
    pt8 = (x_coord + 0.15 * width, y_coord + height)
    points = [pt1, pt2, pt3, pt4, pt5, pt6, pt7, pt8]
    draw_polygon(output, points, colour, opacity)
    if height > 8:
        pt1 = (x_coord + 0.5 * width, y_coord + 0.2 * height)
        pt2 = (x_coord + 0.34 * width, y_coord + 0.6 * height - 1)
        pt3 = (x_coord + 0.64 * width, y_coord + 0.6 * height - 1)
        draw_polygon(output, [pt1, pt2, pt3], WHITE)
    output.write("<!-- End 'A' -->\n")


def draw_letter_c(output, x_coord, y_coord, width, height, colour=BLUE,
                  opacity=1.0):
    """
    Print the svg instructions to draw the letter 'C' to the output.

    :arg output: Stream where to write.
    :type output: file
    :arg x_coord: x-coordinate of the bottom left corner of the box containing
        the 'C'.
    :type x_coord: float
    :arg y_coord: y-coordinate of the bottom left corner ot the box containing
        the 'C'.
    :type y_coord: float
    :arg width: Width of the box containing the 'C'.
    :type width: float
    :arg height: Height of the box containing the 'C'.
    :type height: float
    :arg colour: Color of the box containing the 'C' (default: black). The
        colour is given in RGB.
    :type colour: tuple (int, int, int)
    :arg opacity: Opacity of the box containing the 'C' (default: 1.0)
    :type opacity: float

    :note: The 'C' is drawn by creating the surrounding ellipse with the right
       colour and then superposing an smaller white ellipse and a white
       rectangle on the right side of the ellipe to create the 'C'.

    """

    output.write("<!-- Begin 'C' -->\n")
    draw_ellipse(output, x_coord + width / 2, y_coord + height / 2, width / 2,
                 height / 2, colour, opacity)
    draw_ellipse(output, x_coord + width / 2, y_coord + height / 2, width / 3,
                 height / 3, WHITE)
    draw_rectangle(output, x_coord + width / 2, y_coord + height / 4,
                   0.5 * height, width / 2, WHITE)
    output.write("<!-- End 'C' -->\n")


def draw_letter_g(output, x_coord, y_coord, width, height, colour=YELLOW,
                  opacity=1.0):
    """
    Print the svg instructions to draw the letter 'G' to the output.

    :arg output: Stream where to write.
    :type output: file
    :arg x_coord: x-coordinate of the bottom left corner of the box containing
        the 'G'.
    :type x_coord: float
    :arg y_coord: y-coordinate of the bottom left corner ot the box containing
        the 'G'.
    :type y_coord: float
    :arg width: Width of the box containing the 'G'.
    :type width: float
    :arg height: Height of the box containing the 'G'.
    :type height: float
    :arg colour: Color of the box containing the 'G' (default: black). The
        colour is given in RGB.
    :type colour: tuple (int, int, int)
    :arg opacity: Opacity of the box containing the 'G' (default: 1.0)
    :type opacity: float

    :note: The 'G' is drawn by creating a 'C' and then adding the two other
       rectangles.

    """

    output.write("<!-- Begin 'G' -->\n")
    draw_ellipse(output, x_coord + width / 2, y_coord + height / 2, width / 2,
                 height / 2, colour, opacity)
    draw_ellipse(output, x_coord + width / 2, y_coord + height / 2, width / 3,
                 height / 3, WHITE)
    draw_rectangle(output, x_coord + width / 2, y_coord + height / 4, height /
                   2, width / 2, WHITE)
    draw_rectangle(output, x_coord + width / 2, y_coord + 5 * height / 8,
                   height / 8, width / 2, colour, opacity)
    draw_rectangle(output, x_coord + 7 * width / 8, y_coord + 5 * height / 8,
                   height - 5 * height / 8, width / 8, colour, opacity)
    output.write("<!-- End 'G' -->\n")


def draw_letter_t(output, x_coord, y_coord, width, height, colour=GREEN,
                  opacity=1.0):
    """
    Print the svg instructions to draw the letter 'T' to the output.

    :arg output: Stream where to write.
    :type output: file
    :arg x_coord: x-coordinate of the bottom left corner of the box containing
        the 'T'.
    :type x_coord: float
    :arg y_coord: y-coordinate of the bottom left corner ot the box containing
        the 'T'.
    :type y_coord: float
    :arg width: Width of the box containing the 'T'.
    :type width: float
    :arg height: Height of the box containing the 'T'.
    :type height: float
    :arg colour: Color of the box containing the 'T' (default: black). The
        colour is given in RGB.
    :type colour: tuple (int, int, int)
    :arg opacity: Opacity of the box containing the 'T' (default: 1.0)
    :type opacity: float

    :note: The 'T' is drawn by using two rectangles.

    """

    output.write("<!-- Begin 'T' -->\n")
    draw_rectangle(output, x_coord, y_coord, 0.16 * height, width, colour,
                   opacity)
    draw_rectangle(output, x_coord + 0.4 * width, y_coord, height, .2 * width,
                   colour, opacity)
    output.write("<!-- End 'T' -->\n")


def draw_left_letters(output):
    """
    Print the big four nucleotides on the left side of the svg to the output.

    :arg output: Stream where to write.
    :type output: file

    """
    draw_letter_a(output, 10., 10., 60., 100.)
    draw_letter_c(output, 10., 120., 60., 100.)
    draw_letter_g(output, 10., 230., 60., 100.)
    draw_letter_t(output, 10., 340., 60., 100.)


def draw_dense_letters(output, probabilities, width, xposition, yposition,
                       step, previous_position_proba, index_letter):
    """
    Print the svg instructions to draw the four nucleotides in the dense logo
    with respect to their probability of appearance at the current position,
    and return the y-coordinate position for the next box where to draw
    letters.

    :arg output: Stream where to write.
    :type output: file
    :arg probabilities: Probabilities of appearance of the four nucleotides.
        The given probabilities are given for A, C, G, and T in this order.
    :type probabilities: list of tuple (float, str)
    :arg width: Width of box containing the letters to draw.
    :type width: float
    :arg xposition: x-coordinate of the bottom left corner of the box
        containing the letters to draw.
    :type xposition: float
    :arg yposition: y-coordinate of the bottom left corner of the box
        containing the letters to draw.
    :type yposition: float
    :arg step: Distance between two boxes containing the letters for the same
       position (it is a dense logo!).
    :type step: float
    :arg previous_position_proba: Probabilities of getting ACGT at the
        previous position of the motif.
    :type previous_position_proba: dict
    :arg index_letter: Index of the letter in the alphabet from A:0, C:1, G:2,
        T:3 (i.e. -1<index_letter<4).
    :type index_letter: int

    :returns: The y-coordinate for the next box where to draw letters.
    :rtype: float

    """

    for j in xrange(4):
        proba, letter = probabilities[j]
        height = 100. * proba
        if letter == 'A':
            draw_letter_a(
                output, xposition, yposition, width, height,
                opacity=previous_position_proba[index_letter])
        elif letter == 'C':
            draw_letter_c(
                output, xposition, yposition, width, height,
                opacity=previous_position_proba[index_letter])
        elif letter == 'G':
            draw_letter_g(
                output, xposition, yposition, width, height,
                opacity=previous_position_proba[index_letter])
        elif letter == 'T':
            draw_letter_t(
                output, xposition, yposition, width, height,
                opacity=previous_position_proba[index_letter])
        if j < 3:
            yposition += height
    yposition += height + step
    return yposition


def draw_summary_letters(output, probabilities, width, xposition, yposition,
                         step, entropy):
    """
    Print the svg instruction to writh the ACGT letters at the current position
    and return the x-coordinate and y-coordinate for the box containing the
    letters at the next position.

    :arg output: Stream where to write.
    :type output: file
    :arg probabilities: Probabilities of appearance of the four nucleotides.
        The given probabilities are given for A, C, G, and T in this order.
    :type probabilities: list of tuple (float, str)
    :arg width: Width of box containing the letters to draw.
    :type width: float
    :arg xposition: x-coordinate of the bottom left corner of the box
        containing the letters to draw.
    :type xposition: float
    :arg yposition: y-coordinate of the bottom left corner of the box
        containing the letters to draw.
    :type yposition: float
    :arg step: Distance between two boxes containing the letters for the same
       position (it is a dense logo!).
    :type step: float
    :arg entropy: Entropy of the current position to determine the height of
       the letters.
    :type entropy: float

    :returns: The coordinates (x, y) of the box that will be drawn with letters
        at the next position.
    :rtype: tuple (float, float).

    """

    max_size = (2. - entropy) * 50.
    yposition = 10. + 100. - max_size
    for j in range(0, 4):
        proba, letter = probabilities[j]
        height = max_size * proba
        if letter == 'A':
            draw_letter_a(output, xposition, yposition, width, height)
        elif letter == 'C':
            draw_letter_c(output, xposition, yposition, width, height)
        elif letter == 'G':
            draw_letter_g(output, xposition, yposition, width, height)
        elif letter == 'T':
            draw_letter_t(output, xposition, yposition, width, height)
        if j < 3:
            yposition += height
    return xposition + width + step, yposition + height + step


def draw_frame(output, header_width, header_height, axis_length, axis_width,
               axis_height):
    """
    Print the svg instructions to draw the frame of the logo (i.e. Y- and
    X-axis).

    :arg output: Stream where to write.
    :type output: file
    :arg header_width: Width of the svg.
    :type header_width: float
    :arg header_height: Height of the svg.
    :type header_height: float
    :arg axis_length: number of positions in the motif.
    :type axis_length: float
    :arg axis_width: Width of the logo x-axis.
    :type axis_width: float
    :arg axis_height: Height of the logo y-axis.
    :type axis_height: float

    """

    print_header(output, header_width, header_height)
    draw_axis(output, axis_length, axis_width, axis_height)


def draw_logo_letters(output, tffm, logo_type, xposition=90., width=50.,
                      step=10.):
    """
    Print the svg instructions to draw the letters of the logo and return the
    information content to be printed.

    :arg output: Stream where to write the svg instruction to draw the logo.
    :type output: file
    :arg tffm: The TFFM for which drawing the logo.
    :type tffm: :class:`TFFM`
    :arg logo_type: Kind of logo to draw (either 'summary' or 'dense')
    :type logo_type: str
    :arg xposition: x-coordinate where to start the logo (default: 90.)
    :type xposition: float
    :arg width: Width of the logo.
    :type within: float
    :arg step: Distance between two boxes containing letters in the logo.
    :type step: float

    :note: The computation of the information content is done within the
        drawing since it follows the same algorithm computing the emission
        probabilities. So we do not call the get_information_content method for
        an algorithmic improvement.

    """
    previous_position_proba = tffm.background_emission_proba()
    information_content = 0.
    start = tffm.get_position_start()
    for position in range(start, len(tffm) + start):
        position_proba = {'A': 0., 'C': 0., 'G': 0., 'T': 0.}
        yposition = 10.
        if tffm.kind == TFFM_KIND.ZERO_ORDER:
            __ = tffm.get_emission_update_pos_proba(position_proba, position,
                                                    previous_position_proba, 0)
        else:
            for i in range(0, 4):
                if logo_type == LOGO_TYPE.SUMMARY:
                    __ = tffm.get_emission_update_pos_proba(position_proba,
                                                            position,
                                                            previous_position_proba,
                                                            i)
                else:
                    emissions = tffm.get_emission_update_pos_proba(
                        position_proba, position, previous_position_proba, i)
                    yposition = draw_dense_letters(output, emissions, width,
                                                   xposition, yposition, step,
                                                   previous_position_proba,
                                                   ALPHABET[i])
        previous_position_proba = position_proba.copy()
        if tffm.kind == TFFM_KIND.DETAILED:
            somme = sum(previous_position_proba.values())
            for nucleotide in iter(previous_position_proba):
                previous_position_proba[nucleotide] /= somme
        values = previous_position_proba.values()
        keys = previous_position_proba.keys()
        emissions_tuple = zip(values, keys)
        entropy = tffm_module.compute_entropy(emissions_tuple)
        information_content += 2. - entropy
        if logo_type == LOGO_TYPE.SUMMARY:
            emissions_tuple.sort(reverse=True)
            xposition, yposition = draw_summary_letters(output,
                                                        emissions_tuple, width,
                                                        xposition, yposition,
                                                        step, entropy)
        else:
            xposition += width + step
    return information_content


# TODO raise an error if the logo type is wrong
def draw_logo(output, tffm, logo_type):
    """
    Print the svg instructions to draw the logo to the output.

    :arg output: Stream where to write the svg instruction to draw the logo.
    :type output: file
    :arg tffm: TFFM from which to draw the logo.
    :type tffm: :class:`TFFM`
    :arg logo_type: Kind of logo to draw (either LOGO_TYPE.SUMMARY or
        LOGO_TYPE.DENSE)
    :type logo_type: Enum

    :todo: raise an error if the logo type is wrong.

    """

    header_y_coord = 60 * (len(tffm) + 2)
    if logo_type == LOGO_TYPE.SUMMARY:
        draw_frame(output, 170, header_y_coord, len(tffm), 80, 120)
        draw_y_axis(output)
    else:
        draw_frame(output, 500, header_y_coord, len(tffm), 80, 450)
        draw_left_letters(output)
    information_content = draw_logo_letters(output, tffm, logo_type)
    if logo_type == LOGO_TYPE.SUMMARY:
        print_ic(output, len(tffm) + 2, information_content, 80, 120)
    else:
        print_ic(output, len(tffm) + 2, information_content, 80, 450)
    print_footer(output)
