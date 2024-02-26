#ifndef NCURSES_PLUSPLUS
#define NCURSES_PLUSPLUS

#include <ncurses.h>
#include <vector>
#include <utility>
#include <locale.h>
#include <string>
#include <math.h>
#include <cmath>
#include <ctime>

/// @brief The color pair for black on black (assuming unchanged default color initiations)
#define NPP_BLACK 0
/// @brief The color pair for white on black (assuming unchanged default color initiations)
#define NPP_WHITE 1
/// @brief The color pair for light gray on black (assuming unchanged default color initiations)
#define NPP_LIGHT_GRAY 2
/// @brief The color pair for light grey on black (assuming unchanged default color initiations)
#define NPP_LIGHT_GREY 2
/// @brief The color pair for dark gray on black (assuming unchanged default color initiations)
#define NPP_DARK_GRAY 3
/// @brief The color pair for dark grey on black (assuming unchanged default color initiations)
#define NPP_DARK_GREY 3
/// @brief The color pair for brown on black (assuming unchanged default color initiations)
#define NPP_BROWN 4
/// @brief The color pair for red on black (assuming unchanged default color initiations)
#define NPP_RED 5
/// @brief The color pair for maroon on black (assuming unchanged default color initiations)
#define NPP_MAROON 6
/// @brief The color pair for orange on black (assuming unchanged default color initiations)
#define NPP_ORANGE 7
/// @brief The color pair for yellow on black (assuming unchanged default color initiations)
#define NPP_YELLOW 8
/// @brief The color pair for lime on black (assuming unchanged default color initiations)
#define NPP_LIME 9
/// @brief The color pair for green on black (assuming unchanged default color initiations)
#define NPP_GREEN 10
/// @brief The color pair for cyan on black (assuming unchanged default color initiations)
#define NPP_CYAN 11
/// @brief The color pair for teal on black (assuming unchanged default color initiations)
#define NPP_TEAL 12
/// @brief The color pair for blue on black (assuming unchanged default color initiations)
#define NPP_BLUE 13
/// @brief The color pair for magenta on black (assuming unchanged default color initiations)
#define NPP_MAGENTA 14
/// @brief The color pair for purple on black (assuming unchanged default color initiations)
#define NPP_PURPLE 15

/// @brief Unknown mouse input
#define M_UNKNOWN -1
/// @brief Mouse 1 (Left Click) being pressed
#define M1_PRESS 0
/// @brief Mouse 1 (Left Click) being released
#define M1_RELEASE 1
/// @brief Mouse 1 (Left Click) being clicked once quickly
#define M1_CLICK 2
/// @brief Mouse 1 (Left Click) being clicked twice quickly
#define M1_DCLICK 3
/// @brief Mouse 1 (Left Click) being clicked thrice quickly
#define M1_TCLICK 4
/// @brief Mouse 2 (Right Click) being pressed
#define M2_PRESS 5
/// @brief Mouse 2 (Right Click) being released
#define M2_RELEASE 6
/// @brief Mouse 2 (Right Click) being clicked once quickly
#define M2_CLICK 7
/// @brief Mouse 2 (Right Click) being clicked twice quickly
#define M2_DCLICK 8
/// @brief Mouse 2 (Right Click) being clicked thrice quickly
#define M2_TCLICK 9
/// @brief Mouse 3 (Middle Click) being pressed
#define M3_PRESS 10
/// @brief Mouse 3 (Middle Click) being released
#define M3_RELEASE 11
/// @brief Mouse 3 (Middle Click) being clicked once quickly
#define M3_CLICK 12
/// @brief Mouse 3 (Middle Click) being clicked twice quickly
#define M3_DCLICK 13
/// @brief Mouse 3 (Middle Click) being clicked thrice quickly
#define M3_TCLICK 14
/// @brief Mouse 4 (Scroll Up) being pressed (ACTUAL INPUT - USE THIS FOR SCROLLING)
#define M4_PRESS 15
/// @brief Mouse 4 (Scroll Up) being released (DOESN'T DO ANYTHING)
#define M4_RELEASE 16
/// @brief Mouse 4 (Scroll Up) being clicked once quickly (DOESN'T DO ANYTHING)
#define M4_CLICK 17
/// @brief Mouse 4 (Scroll Up) being clicked twice quickly (DOESN'T DO ANYTHING)
#define M4_DCLICK 18
/// @brief Mouse 4 (Scroll Up) being clicked thrice quickly (DOESN'T DO ANYTHING)
#define M4_TCLICK 19
/// @brief Mouse 5 (Scroll Down) being pressed (ACTUAL INPUT - USE THIS FOR SCROLLING)
#define M5_PRESS 20
/// @brief Mouse 5 (Scroll Down) being released (DOESN'T DO ANYTHING)
#define M5_RELEASE 21
/// @brief Mouse 5 (Scroll Down) being clicked once quickly (DOESN'T DO ANYTHING)
#define M5_CLICK 22
/// @brief Mouse 5 (Scroll Down) being clicked twice quickly (DOESN'T DO ANYTHING)
#define M5_DCLICK 23
/// @brief Mouse 5 (Scroll Down being clicked thrice quickly (DOESN'T DO ANYTHING)
#define M5_TCLICK 24

/// @brief Light box drawing characters with square corners
#define LIGHT_HARD 0
/// @brief Light box drawing characters with rounded corners
#define LIGHT_SOFT 1
/// @brief Heavy box drawing characters in both the horiztonal and vertical directions
#define HEAVY_BOTH 2
/// @brief Heavy box drawing characters in the horiztonal direction, light in the vertical directions
#define HEAVY_HORIZONTAL 3
/// @brief Heavy box drawing characters in the vertical direction, light in the horizontal directions
#define HEAVY_VERTICAL 4
/// @brief Doubled box drawing characters in both the horizontal and vertical directions
#define DOUBLED_BOTH 5
/// @brief Doubled box drawing characters in the horiztonal direction, light in the vertical directions
#define DOUBLED_HORIZTONAL 6
/// @brief Doubled box drawing characters in the vertical direction, light in the horizontal directions
#define DOUBLED_VERTICAL 7

/// @brief Solid Line
#define DASHED_NONE 0
/// @brief Dashed line with one gap per character
#define DASHED_DOUBLE 1
/// @brief Dashed line with two gaps per character
#define DASHED_TRIPLE 2
/// @brief Dashed line with three gap per character
#define DASHED_QUADRUPLE 3
/// @brief Half of a solid line used as a dashed line (left/up)
#define DASHED_SINGLE_1 4
/// @brief Half of a solid line used as a dashed line (right/down)
#define DASHED_SINGLE_2 5

/// @brief Matrix Text in a 6x6 grid (3x3 characters)
#define MTEXT_6x6 0
/// @brief Matrix Text in an 8x8 grid (4x4 characters)
#define MTEXT_8x8 1

/// @brief Any character that isn't a box-drawing one
#define STYLE_NONE 0
/// @brief Box-drawing characters with a light line weight
#define STYLE_LIGHT 1
/// @brief Box-drawing characters with a heavy line weight
#define STYLE_HEAVY 2
/// @brief Box-drawing characters with two lightweight lines (doubled lines)
#define STYLE_DOUBLED 3

/// @brief The top edge of box-drawing characters
#define DIR_UP 0
/// @brief The bottom edge of box-drawing characters
#define DIR_DOWN 1
/// @brief The left edge of box drawing characters
#define DIR_LEFT 2
/// @brief The right edge of box-drawing characters
#define DIR_RIGHT 3

/// @brief ncursespp - A namespace that provides alternate classes/functions to the ones from ncurses.h designed to make programming TUIs easier; there is a decent chance that some functionality is missing from the original header file
namespace npp {
    /// @brief Default values that can be changed in place of providing parameter arguments for many functions
    static struct {
        /// @brief Color pair
        unsigned char Color = NPP_WHITE;
        /// @brief The set of attributes to use by default when writing to the window
        std::string Attributes = "";
        /// @brief The pair of y (row) and x (col) offsets to use when chaining writing functions
        std::pair<unsigned short, unsigned short> Offset = {0, 1};
        /// @brief Whether to wrap to the starting x or to the window's x
        bool Wrap = true;
        /// @brief The minimum width of various numbers to have (in characters) when written to the window
        unsigned short MinWidth = 0;

        unsigned char Matrix = 0;
        bool DanglingLetters = true;
        unsigned short Kerning = 1;
        unsigned short Leading = 1;

        bool UsePadding = true;

        std::pair<unsigned char, unsigned char> Style = {LIGHT_HARD, DASHED_NONE};
        bool Mergeable = true;
        bool CanMerge = true;
    } Defaults;

    /// @brief Various keys used to store and request information relevant to each key
    static const struct {
        /// @brief The key that matches unicode box drawing characters to the directions/styles that make them up
        std::vector<std::pair<wchar_t, std::vector<unsigned char>>> LineDrawing = {
            // Horizontal and vertical lines off all 6 dash styles (none, single 1, single 2, double, triple, quadruple)
            {L'─', {0, 0, 1, 1}}, {L'╌', {0, 0, 1, 1}}, {L'┄', {0, 0, 1, 1}}, {L'┈', {0, 0, 1, 1}}, {L'╴', {0, 0, 1, 1}}, {L'╶', {0, 0, 1, 1}},
            {L'━', {0, 0, 2, 2}}, {L'╍', {0, 0, 2, 2}}, {L'┅', {0, 0, 2, 2}}, {L'┉', {0, 0, 2, 2}}, {L'╸', {0, 0, 2, 2}}, {L'╺', {0, 0, 2, 2}},
            {L'│', {1, 1, 0, 0}}, {L'╎', {1, 1, 0, 0}}, {L'┆', {1, 1, 0, 0}}, {L'┊', {1, 1, 0, 0}}, {L'╵', {1, 1, 0, 0}}, {L'╷', {1, 1, 0, 0}},
            {L'┃', {2, 2, 0, 0}}, {L'╏', {2, 2, 0, 0}}, {L'┇', {2, 2, 0, 0}}, {L'┋', {2, 2, 0, 0}}, {L'╹', {2, 2, 0, 0}}, {L'╻', {2, 2, 0, 0}},
            {L'═', {0, 0, 3, 3}}, {L'║', {3, 3, 0, 0}},
            // Corners
            {L'┌', {0, 1, 0, 1}}, {L'╭', {0, 1, 0, 1}}, {L'┐', {0, 1, 1, 0}}, {L'╮', {0, 1, 1, 0}}, {L'└', {1, 0, 0, 1}}, {L'╰', {1, 0, 0, 1}}, {L'┘', {1, 0, 1, 0}}, {L'╯', {1, 0, 1, 0}}, 
            {L'┍', {0, 1, 0, 2}}, {L'┎', {0, 2, 0, 1}}, {L'┏', {0, 2, 0, 2}}, {L'╒', {0, 1, 0, 3}}, {L'╓', {0, 3, 0, 1}}, {L'╔', {0, 3, 0, 3}},
            {L'┑', {0, 1, 2, 0}}, {L'┒', {0, 2, 1, 0}}, {L'┓', {0, 2, 2, 0}}, {L'╕', {0, 1, 3, 0}}, {L'╖', {0, 3, 1, 0}}, {L'╗', {0, 3, 3, 0}},
            {L'┕', {1, 0, 0, 2}}, {L'┖', {2, 0, 0, 1}}, {L'┗', {2, 0, 0, 2}}, {L'╘', {1, 0, 0, 3}}, {L'╙', {3, 0, 0, 1}}, {L'╚', {3, 0, 0, 3}},
            {L'┙', {1, 0, 2, 0}}, {L'┚', {2, 0, 1, 0}}, {L'┛', {2, 0, 2, 0}}, {L'╛', {1, 0, 3, 0}}, {L'╜', {3, 0, 1, 0}}, {L'╝', {3, 0, 3, 0}},
            // Tees
            {L'├', {1, 1, 0, 1}}, {L'┝', {1, 1, 0, 2}}, {L'┞', {2, 1, 0, 1}}, {L'┟', {1, 2, 0, 1}}, {L'┠', {2, 2, 0, 1}}, {L'┡', {2, 1, 0, 2}}, {L'┢', {1, 2, 0, 2}}, {L'┣', {2, 2, 0, 2}}, {L'╞', {1, 1, 0, 3}}, {L'╟', {3, 3, 0, 1}}, {L'╠', {3, 3, 0, 3}},
            {L'┤', {1, 1, 1, 0}}, {L'┥', {1, 1, 2, 0}}, {L'┦', {2, 1, 1, 0}}, {L'┧', {1, 2, 1, 0}}, {L'┨', {2, 2, 1, 0}}, {L'┩', {2, 1, 2, 0}}, {L'┪', {1, 2, 2, 0}}, {L'┫', {2, 2, 2, 0}}, {L'╡', {1, 1, 3, 0}}, {L'╢', {3, 3, 1, 0}}, {L'╣', {3, 3, 3, 0}},
            {L'┬', {0, 1, 1, 1}}, {L'┭', {0, 1, 2, 1}}, {L'┮', {0, 1, 1, 2}}, {L'┯', {0, 1, 2, 2}}, {L'┰', {0, 2, 1, 1}}, {L'┱', {0, 2, 2, 1}}, {L'┲', {0, 2, 1, 2}}, {L'┳', {0, 2, 2, 2}}, {L'╤', {0, 1, 3, 3}}, {L'╥', {0, 3, 1, 1}}, {L'╦', {0, 3, 3, 3}},
            {L'┴', {1, 0, 1, 1}}, {L'┵', {1, 0, 2, 1}}, {L'┶', {1, 0, 1, 2}}, {L'┷', {1, 0, 2, 2}}, {L'┸', {2, 0, 1, 1}}, {L'┹', {2, 0, 2, 1}}, {L'┺', {2, 0, 1, 2}}, {L'┻', {2, 0, 2, 2}}, {L'╧', {1, 0, 3, 3}}, {L'╨', {3, 0, 1, 1}}, {L'╩', {3, 0, 3, 3}},
            // Crosses
            {L'┼', {1, 1, 1, 1}}, {L'┽', {1, 1, 2, 1}}, {L'┾', {1, 1, 1, 2}}, {L'┿', {1, 1, 2, 2}}, {L'╀', {2, 1, 1, 1}}, {L'╁', {1, 2, 1, 1}}, {L'╂', {2, 2, 1, 1}}, {L'╃', {2, 1, 2, 1}},
            {L'╄', {2, 1, 1, 2}}, {L'╅', {1, 2, 2, 1}}, {L'╆', {1, 2, 1, 2}}, {L'╇', {2, 2, 1, 2}}, {L'╈', {1, 2, 2, 2}}, {L'╉', {2, 2, 2, 1}}, {L'╊', {2, 2, 1, 2}}, {L'╋', {2, 2, 2, 2}},
            {L'╪', {1, 1, 3, 3}}, {L'╫', {3, 3, 1, 1}}, {L'╬', {3, 3, 3, 3}},
            // Thickness transitions
            {L'╼', {0, 0, 1, 2}}, {L'╽', {1, 2, 0, 0}}, {L'╾', {0, 0, 2, 1}}, {L'╿', {2, 1, 0, 0}}
        };

        /// @brief The key that matches a 3x3 set of unicode box element characters to the regular character they represent
        std::vector<std::pair<char, std::vector<std::vector<std::wstring>>>> MatrixText = {
            {' ', {{L"   ", L"   ", L"   "}, {L"    ", L"    ", L"    ", L"    "}}},
            {'a', {{L"   ", L"▞▀▟", L"▚▄▜"}, {L"▗▄▄▖", L"▗▄▄▐", L"▌  █", L"▚▄▄▜"}}},
            {'b', {{L"▌  ", L"▙▀▚", L"▛▄▞"}, {L"▌   ", L"▌▄▄▖", L"█  ▐", L"▛▄▄▞"}}},
            {'c', {{L"   ", L"▞▀▀", L"▚▄▄"}, {L"    ", L"▗▄▄▄", L"▌   ", L"▚▄▄▄"}}},
            {'d', {{L"  ▐", L"▞▀▟", L"▚▄▜"}, {L"   ▐", L"▗▄▄▐", L"▌  █", L"▚▄▄▜"}}},
            {'e', {{L"   ", L"▟█▙", L"▚▄▄"}, {L"    ", L"▗▄▄▖", L"▙▄▄▟", L"▚▄▄▄"}}},
            {'f', {{L" ▞▖", L"▗▙▖", L" ▌ "}, {L"  ▞▖", L" ▄▙▖", L"  ▌ ", L"  ▚ "}}},
            {'g', {{L"▞▀▟", L"▚▄▜", L"▗▄▞"}, {L"▗▄▄▗", L"▌  █", L"▚▄▄▜", L"▗▄▄▞"}}},
            {'h', {{L"▌  ", L"▙▀▚", L"▌ ▐"}, {L"▌   ", L"▌▄▄▖", L"█  ▐", L"▌  ▐"}}},
            {'i', {{L" ▘ ", L" ▌ ", L" ▚ "}, {L"  ▖ ", L"  ▖ ", L"  ▌ ", L"  ▚ "}}},
            {'j', {{L" ▝ ", L" ▐ ", L"▝▞ "}, {L"  ▖ ", L"  ▖ ", L"  ▌ ", L" ▚▘ "}}},
            {'k', {{L"▌  ", L"▙▄▘", L"▌ ▌"}, {L"▌   ", L"▌  ▗", L"▙▄▄▘", L"▌  ▚"}}},
            {'l', {{L" ▌ ", L" ▌ ", L" ▚ "}, {L"  ▌ ", L"  ▌ ", L"  ▌ ", L"  ▚ "}}},
            {'m', {{L"   ", L"▛▞▖", L"▌▌▌"}, {L"    ", L"▖▄▗▖", L"▛ ▌▐", L"▌ ▌▐"}}},
            {'n', {{L"   ", L"▛▀▚", L"▌ ▐"}, {L"    ", L"▖▄▄▖", L"▛  ▐", L"▌  ▐"}}},
            {'o', {{L"   ", L"▞▀▚", L"▚▄▞"}, {L"    ", L"▗▄▄▖", L"▌  ▐", L"▚▄▄▞"}}},
            {'p', {{L"▞▀▚", L"▙▄▞", L"▌  "}, {L"▖▄▄▖", L"█  ▐", L"▛▄▄▞", L"▌   "}}},
            {'q', {{L"▞▀▚", L"▚▄▟", L"  ▐"}, {L"▗▄▄▗", L"▌  █", L"▚▄▄▜", L"   ▐"}}},
            {'r', {{L"   ", L"▙▀▚", L"▌  "}, {L"    ", L"▖▄▄▖", L"▛  ▝", L"▌   "}}},
            {'s', {{L"▗▄▖", L"▚▄▖", L"▗▄▞"}, {L"    ", L"▗▄▄▖", L"▚▄▄▖", L"▗▄▄▞"}}},
            {'t', {{L" ▌ ", L"▀▛▀", L" ▚ "}, {L"  ▌ ", L" ▄▙▖", L"  ▌ ", L"  ▚ "}}},
            {'u', {{L"   ", L"▌ ▐", L"▚▄▟"}, {L"    ", L"▖  ▗", L"▌  ▐", L"▝▄▄▜"}}},
            {'v', {{L"   ", L"▌ ▐", L"▝▄▘"}, {L"    ", L"▖  ▗", L"▚  ▞", L" ▚▞ "}}},
            {'w', {{L"   ", L"▐▐▐", L"▝▞▟"}, {L"    ", L"▖▗ ▗", L"▌▐ ▐", L"▚▞▄▜"}}},
            {'x', {{L"   ", L"▝▄▘", L"▗▀▖"}, {L"    ", L"▗  ▖", L" ▚▞ ", L"▗▘▝▖"}}},
            {'y', {{L"▌ ▐", L"▚▄▟", L"▗▄▞"}, {L"▖  ▗", L"▌  ▐", L"▝▄▄▌", L"▗▄▄▘"}}},
            {'z', {{L"▄▄▄", L"▗▄▞", L"▙▄▄"}, {L"    ", L"▄▄▄▄", L" ▄▄▘", L"▟▄▄▄"}}},
            {'A', {{L"▞▀▚", L"▙▄▟", L"▌ ▐"}, {L"▞▀▀▚", L"▌  ▐", L"▛▀▀▜", L"▌  ▐"}}},
            {'B', {{L"▛▀▚", L"▛▀▚", L"▙▄▞"}, {L"▛▀▀▚", L"▙▄▄▞", L"▌  ▐", L"▙▄▄▞"}}},
            {'C', {{L"▞▀▚", L"▌  ", L"▚▄▞"}, {L"▞▀▀▚", L"▌   ", L"▌   ", L"▚▄▄▞"}}},
            {'D', {{L"▛▀▚", L"▌ ▐", L"▙▄▞"}, {L"▛▀▀▚", L"▌  ▐", L"▌  ▐", L"▙▄▄▞"}}},
            {'E', {{L"▛▀▀", L"▛▀▀", L"▙▄▄"}, {L"▛▀▀▀", L"▙▄▄▄", L"▌   ", L"▙▄▄▄"}}},
            {'F', {{L"▛▀▀", L"▛▀▀", L"▌  "}, {L"▛▀▀▀", L"▙▄▄▄", L"▌   ", L"▌   "}}},
            {'G', {{L"▞▀▚", L"▌ ▄", L"▚▄▟"}, {L"▞▀▀▚", L"▌   ", L"▌ ▀▜", L"▚▄▄▜"}}},
            {'H', {{L"▌ ▐", L"▛▀▜", L"▌ ▐"}, {L"▌  ▐", L"▙▄▄▟", L"▌  ▐", L"▌  ▐"}}},
            {'I', {{L"▀▜▀", L" ▐ ", L"▄▟▄"}, {L"▀▀▛▀", L"  ▌ ", L"  ▌ ", L"▄▄▙▄"}}},
            {'J', {{L"▀▜▀", L" ▐ ", L"▚▟ "}, {L"▀▀▛▀", L"  ▌ ", L"  ▌ ", L"▚▄▘ "}}},
            {'K', {{L"▌ ▞", L"▛▀▖", L"▌ ▐"}, {L"▌  ▐", L"▙▄▞▘", L"▌ ▝▚", L"▌  ▐"}}},
            {'L', {{L"▌  ", L"▌  ", L"▙▄▄"}, {L"▌   ", L"▌   ", L"▌   ", L"▙▄▄▄"}}},
            {'M', {{L"▙ ▟", L"▌▀▐", L"▌ ▐"}, {L"▙  ▟", L"▌▚▞▐", L"▌  ▐", L"▌  ▐"}}},
            {'N', {{L"▙ ▐", L"▌▚▐", L"▌ ▜"}, {L"▙  ▐", L"▌▚ ▐", L"▌ ▚▐", L"▌  ▜"}}},
            {'O', {{L"▞▀▚", L"▌ ▐", L"▚▄▞"}, {L"▞▀▀▚", L"▌  ▐", L"▌  ▐", L"▚▄▄▞"}}},
            {'P', {{L"▛▀▚", L"▙▄▞", L"▌  "}, {L"▛▀▀▚", L"▙▄▄▞", L"▌   ", L"▌   "}}},
            {'Q', {{L"▞▀▚", L"▌▗▐", L"▚▄▚"}, {L"▞▀▀▚", L"▌  ▐", L"▌ ▚▐", L"▚▄▄▚"}}},
            {'R', {{L"▛▀▚", L"▙▄▞", L"▌ ▐"}, {L"▛▀▀▚", L"▙▄▄▞", L"▌  ▚", L"▌  ▐"}}},
            {'S', {{L"▞▀▘", L"▝▀▚", L"▚▄▞"}, {L"▞▀▀▚", L"▚▄  ", L"  ▀▚", L"▚▄▄▞"}}},
            {'T', {{L"▀▜▀", L" ▐ ", L" ▐ "}, {L"▀▀▛▀", L"  ▌ ", L"  ▌ ", L"  ▌ "}}},
            {'U', {{L"▌ ▐", L"▌ ▐", L"▚▄▞"}, {L"▌  ▐", L"▌  ▐", L"▌  ▐", L"▚▄▄▞"}}},
            {'V', {{L"▌ ▐", L"▚ ▞", L"▝▄▘"}, {L"▌  ▐", L"▌  ▐", L"▚  ▞", L" ▚▞ "}}},
            {'W', {{L"▌ ▐", L"▌▄▐", L"▛ ▜"}, {L"▌  ▐", L"▌  ▐", L"▌▞▚▐", L"▛  ▜"}}},
            {'X', {{L"▚ ▞", L" █ ", L"▞ ▚"}, {L"▚  ▞", L" ▚▞ ", L" ▞▚ ", L"▞  ▚"}}},
            {'Y', {{L"▌ ▐", L"▝▄▘", L" █ "}, {L"▌  ▐", L"▝▖▗▘", L" ▝▌ ", L"  ▌ "}}},
            {'Z', {{L"▀▀▜", L"▗▞▘", L"▙▄▄"}, {L"▀▀▀▜", L"  ▄▘", L"▗▀  ", L"▙▄▄▄"}}},
            {'0', {{L"▞▀▙", L"▌▞▐", L"▜▄▞"}, {L"▞▀▀▙", L"▌ ▞▐", L"▌▞ ▐", L"▜▄▄▞"}}},
            {'1', {{L" ▟ ", L" ▐ ", L" ▟▖"}, {L" ▞▌ ", L"  ▌ ", L"  ▌ ", L" ▄▙▖"}}},
            {'2', {{L"▞▀▚", L" ▗▞", L"▟▙▄"}, {L"▞▀▀▚", L"   ▞", L" ▄▀ ", L"▟▄▄▄"}}},
            {'3', {{L"▞▀▚", L" ▀▚", L"▚▄▞"}, {L"▞▀▀▚", L" ▄▄▞", L"   ▐", L"▚▄▄▞"}}},
            {'4', {{L"▌ ▌", L"▝▀▛", L"  ▌"}, {L"▌  ▌", L"▙▄▄▙", L"   ▌", L"   ▌"}}},
            {'5', {{L"▛▀▀", L"▀▀▚", L"▚▄▞"}, {L"▛▀▀▀", L"▚▄▄▖", L"   ▐", L"▚▄▄▞"}}},
            {'6', {{L"▞▀▀", L"▛▀▚", L"▚▄▞"}, {L"▞▀▀▀", L"▙▄▄▖", L"▌  ▐", L"▚▄▄▞"}}},
            {'7', {{L"▀▀▜", L" ▗▘", L" ▌ "}, {L"▀▀▀▜", L" ▄▄▙", L"  ▌ ", L" ▐  "}}},
            {'8', {{L"▞▀▚", L"▞▀▚", L"▚▄▞"}, {L"▞▀▀▚", L"▚▄▄▞", L"▌  ▐", L"▚▄▄▞"}}},
            {'9', {{L"▞▀▚", L"▚▄▟", L"▗▄▟"}, {L"▞▀▀▚", L"▚▄▄▟", L"   ▐", L"▗▄▄▟"}}},
            {'`', {{L"▗  ", L" ▘ ", L"   "}, {L" ▗  ", L"  ▘ ", L"    ", L"    "}}},
            {'~', {{L"▗▖▗", L"▘▝▘", L"   "}, {L" ▄ ▖", L"▝ ▀ ", L"    ", L"    "}}},
            {'!', {{L" █ ", L" █ ", L" ▄ "}, {L" ▐▌ ", L" ▐▌ ", L" ▝▘ ", L" ▐▌ "}}},
            {'@', {{L"▞▀▚", L"▌█▟", L"▚▄▄"}, {L"▞▀▀▚", L"▌▞▚▐", L"▌▚▟▟", L"▚▄▄▄"}}},
            {'#', {{L"▟▄▙", L"▐ ▌", L"▜▀▛"}, {L" ▌▐ ", L"▀▛▜▀", L"▄▙▟▄", L" ▌▐ "}}},
            {'$', {{L"▗▙▖", L"▚▙▖", L"▗▙▞"}, {L"  ▖ ", L"▞▀▛▘", L"▝▀▛▚", L"▝▀▛▘"}}},
            {'%', {{L"█ ▞", L" ▞ ", L"▞ █"}, {L"▞▚ ▞", L"▚▞▞ ", L" ▞▞▚", L"▞ ▚▞"}}},
            {'^', {{L" ▄ ", L"▝ ▘", L"   "}, {L" ▗▖ ", L" ▘▝ ", L"    ", L"    "}}},
            {'&', {{L"▞▚ ", L"▞▌▖", L"▚▞▖"}, {L"▗▀▖ ", L"▝▄▘ ", L"▞▝▖▐", L"▚▄▞▚"}}},
            {'*', {{L"▝▄▘", L"▗▀▖", L"   "}, {L" ▚▙▘", L" ▘▘▘", L"    ", L"    "}}},
            {'(', {{L" ▞▘", L"▐  ", L" ▚▖"}, {L" ▗▀ ", L" ▌  ", L" ▌  ", L" ▝▄ "}}},
            {')', {{L"▝▚ ", L"  ▌", L"▗▞ "}, {L" ▀▖ ", L"  ▐ ", L"  ▐ ", L" ▄▘ "}}},
            {'-', {{L"   ", L"▗▄▖", L"   "}, {L"    ", L"▗▄▄▖", L"    ", L"    "}}},
            {'_', {{L"   ", L"   ", L"▄▄▄"}, {L"    ", L"    ", L"    ", L"▄▄▄▄"}}},
            {'=', {{L"   ", L"▄█▄", L" ▀ "}, {L"    ", L"▗▄▄▖", L"▗▄▄▖", L"    "}}},
            {'+', {{L"   ", L"▝▀▘", L"▝▀▘"}, {L"  ▖ ", L"▗▄▙▖", L"  ▌ ", L"    "}}},
            {'[', {{L"▐▀▘", L"▐  ", L"▐▄▖"}, {L" ▛▀ ", L" ▌  ", L" ▌  ", L" ▙▄ "}}},
            {'{', {{L" ▛▘", L"█  ", L" ▙▖"}, {L" ▛▀ ", L"▗▘  ", L"▝▖  ", L" ▙▄ "}}},
            {']', {{L"▝▀▌", L"  ▌", L"▗▄▌"}, {L" ▀▜ ", L"  ▐ ", L"  ▐ ", L" ▄▟ "}}},
            {'}', {{L"▝▜ ", L"  █", L"▗▟ "}, {L" ▀▜ ", L"  ▝▖", L"  ▗▘", L" ▄▟ "}}},
            {'\\', {{L" █ ", L" █ ", L" █ "}, {L"▚   ", L" ▚  ", L"  ▚ ", L"   ▚"}}},
            {'|', {{L"▚  ", L" ▚ ", L"  ▚"}, {L" ▐▌ ", L" ▐▌ ", L" ▐▌ ", L" ▐▌ "}}},
            {';', {{L"   ", L" ▀ ", L" ▜ "}, {L"    ", L" ▗▖ ", L" ▝▘ ", L" ▝▌ "}}},
            {':', {{L"   ", L" ▀ ", L" ▄ "}, {L"    ", L" ▗▖ ", L" ▝▘ ", L" ▐▌ "}}},
            {'\'', {{L" ▗▖", L"  ▘", L"   "}, {L" ▖▖ ", L" ▘▘ ", L"    ", L"    "}}},
            {'"', {{L"▗ ▖", L"▝ ▘", L"   "}, {L"  ▖ ", L"  ▘ ", L"    ", L"    "}}},
            {',', {{L"   ", L"   ", L" ▜ "}, {L"    ", L"    ", L"    ", L" ▝▌ "}}},
            {'<', {{L" ▗▖", L"▐▌ ", L" ▝▘"}, {L"  ▄▖", L"▗▀  ", L"▝▄  ", L"  ▀▘"}}},
            {'.', {{L"   ", L"   ", L" ▄ "}, {L"    ", L"    ", L"    ", L" ▐▌ "}}},
            {'>', {{L"▗▖ ", L" ▐▌", L"▝▘ "}, {L"▗▄  ", L"  ▀▖", L"  ▄▘", L"▝▀  "}}},
            {'?', {{L"▞▀▚", L" ▄▘", L" ▄ "}, {L"   ▞", L"  ▞ ", L" ▞  ", L"▞   "}}},
            {'/', {{L"  ▞", L" ▞ ", L"▞  "}, {L"▞▀▀▚", L"  ▗▞", L" ▐▌ ", L" ▗▖ "}}}
        };
        std::vector<unsigned char> DanglingOffsets = {1, 1};
    } Keys;

    /// @brief Contain the data (character/color/attributes/etc) for each cell
    struct Cell {
        /// @brief Character contained in the cell
        wchar_t Char = L' ';

        /// @brief Color pair used for the cell
        unsigned char Color = 1;

        /// @brief Whether the cell is bolded or not
        bool Bold = false;
        /// @brief Whether the cell is italicized or not
        bool Italic = false;
        /// @brief Whether the cell is underlined or not
        bool Under = false;
        /// @brief Whether the cell is reversed (color pair) or not
        bool Rev = false;
        /// @brief Whether the cell blinks or not
        bool Blink = false;
        /// @brief Whether the cell is dim or not
        bool Dim = false;
        /// @brief Whether the cell is invisible or not
        bool Invis = false;
        /// @brief Whether the cell should stand out or not
        bool Stand = false;
        /// @brief Whether the cell should be protected or not
        bool Prot = false;
        /// @brief Whether the cell should use the alternate character set or not
        bool Alt = false;

        /// @brief If the cell is a line drawing character, whether it can be merged with or not
        bool CanMerge = false;
    };

    /// @brief Essentially the same as the MEVENT from ncurses, but with friendlier values
    class {
        private:
            /// @brief Last button/input the mouse had
            char Button = -1;
            /// @brief x-position (col) of the mouse cursor during the last event
            unsigned short X = 0;
            /// @brief y-position (row) of the mouse cursor during the last event
            unsigned short Y = 0;
            /// @brief z-position (???) of the mouse cursor during the last event
            unsigned short Z = 0;

        public:
            /// @brief Get Input - Get the mouse's last recorded input
            /// @returns A number corresponding to the mouse's last recorded input
            const char ginput() {return Button;}
            /// @brief Get X - Get the mouse's last recorded x-position (col)
            /// @returns The mouse's last recorded x-position (col)
            const unsigned short gx() {return X;}
            /// @brief Get Y - Get the mouse's last recorded y-position (row)
            /// @returns The mouse's last recorded y-position (row)
            const unsigned short gy() {return Y;}
            /// @brief Get Z - Get the mouse's last recorded z-position (???)
            /// @returns The mouse's last recorded z-position (???)
            const unsigned short gz() {return Z;}

            /// @brief Get Mouse - Take an input from gchar() or something and update the mouse
            /// @param input Integer input (most likely from gchar())
            /// @returns True if the Mouse was successfully updated, false if not
            bool gmouse(int input) {
                MEVENT event;
                if (getmouse(&event) != OK) {return false;}

                // Change the janky ncurses values for the mouse into easier ones
                // Mouse 2 is actually bound to right click as well, not middle click
                switch (event.bstate) {
                    case BUTTON1_PRESSED:
                        Button = M1_PRESS;
                        break;
                    case BUTTON1_RELEASED:
                        Button = M1_RELEASE;
                        break;
                    case BUTTON1_CLICKED:
                        Button = M1_CLICK;
                        break;
                    case BUTTON1_DOUBLE_CLICKED:
                        Button = M1_DCLICK;
                        break;
                    case BUTTON1_TRIPLE_CLICKED:
                        Button = M1_TCLICK;
                        break;
                    case BUTTON2_PRESSED:
                        Button = M3_PRESS;
                        break;
                    case BUTTON2_RELEASED:
                        Button = M3_RELEASE;
                        break;
                    case BUTTON2_CLICKED:
                        Button = M3_CLICK;
                        break;
                    case BUTTON2_DOUBLE_CLICKED:
                        Button = M3_DCLICK;
                        break;
                    case BUTTON2_TRIPLE_CLICKED:
                        Button = M3_TCLICK;
                        break;
                    case BUTTON3_PRESSED:
                        Button = M2_PRESS;
                        break;
                    case BUTTON3_RELEASED:
                        Button = M2_RELEASE;
                        break;
                    case BUTTON3_CLICKED:
                        Button = M2_CLICK;
                        break;
                    case BUTTON3_DOUBLE_CLICKED:
                        Button = M2_DCLICK;
                        break;
                    case BUTTON3_TRIPLE_CLICKED:
                        Button = M2_TCLICK;
                        break;
                    case BUTTON4_PRESSED:
                        Button = M4_PRESS;
                        break;
                    case BUTTON4_RELEASED:
                        Button = M4_RELEASE;
                        break;
                    case BUTTON4_CLICKED:
                        Button = M4_CLICK;
                        break;
                    case BUTTON4_DOUBLE_CLICKED:
                        Button = M4_DCLICK;
                        break;
                    case BUTTON4_TRIPLE_CLICKED:
                        Button = M4_TCLICK;
                        break;
                    case BUTTON5_PRESSED:
                        Button = M5_PRESS;
                        break;
                    case BUTTON5_RELEASED:
                        Button = M5_RELEASE;
                        break;
                    case BUTTON5_CLICKED:
                        Button = M5_CLICK;
                        break;
                    case BUTTON5_DOUBLE_CLICKED:
                        Button = M5_DCLICK;
                        break;
                    case BUTTON5_TRIPLE_CLICKED:
                        Button = M5_TCLICK;
                        break;
                    default:
                        Button = M_UNKNOWN;
                        break;
                }

                X = event.x;
                Y = event.y;
                Z = event.z;

                return true;
            }
    } Mouse;

    /// @brief The npp version of the WINDOW class from ncurses.h - comes with better support for unicode characters, much better line drawing capabilities, flashy rendering animations, and other fun bonuses
    class Window {
        private:
            /// @brief ncurses window that allows everything to be interacted with inside a compatible terminal
            WINDOW *Win;

            /// @brief Y-dimension (rows) of the window
            unsigned short DimY;
            /// @brief X-dimension (cols) of the window
            unsigned short DimX;
            /// @brief Y-position (row) of the top-left corner of the window
            unsigned short PosY;
            /// @brief X-position (col) of the top-left corner of the window
            unsigned short PosX;

            /// @brief Amount of characters to add padding to from the top of the screen
            unsigned short PadUp = 0;
            /// @brief Amount of characters to add padding to from the bottom of the screen
            unsigned short PadDown = 0;
            /// @brief Amount of characters to add padding to from the left of the screen
            unsigned short PadLeft = 0;
            /// @brief Amount of characters to add padding to from the right of the screen
            unsigned short PadRight = 0;

            /// @brief If the user can skip wait() functions with an input
            bool CanSkip = true;

            std::vector<std::vector<npp::Cell>> Grid;

            //
            // INTERFACING WITH NCURSES
            //

            /// @brief Write - Write a cell to the ncurses window
            /// @param y Y-position (row) of the cell to write
            /// @param x X-position (col) of the cell to write
            void write(unsigned short y, unsigned short x) {
                // No printing characters outside of the ncurses window
                if (y < 0 || y >= DimY || x < 0 || x >= DimX) {return;}

                wattron(Win, COLOR_PAIR(Grid[y][x].Color));
                toggleAttributes(Grid[y][x].Bold, Grid[y][x].Italic, Grid[y][x].Under, Grid[y][x].Rev, Grid[y][x].Blink, Grid[y][x].Dim, Grid[y][x].Invis, Grid[y][x].Stand, Grid[y][x].Prot, Grid[y][x].Alt);

                if (Grid[y][x].Char == '%') {mvwprintw(Win, y, x, "%%");}
                else {
                    std::wstring str; str += Grid[y][x].Char;
                    mvwaddwstr(Win, y, x, str.c_str());
                }

                wattroff(Win, COLOR_PAIR(Grid[y][x].Color));
                toggleAttributes(Grid[y][x].Bold, Grid[y][x].Italic, Grid[y][x].Under, Grid[y][x].Rev, Grid[y][x].Blink, Grid[y][x].Dim, Grid[y][x].Invis, Grid[y][x].Stand, Grid[y][x].Prot, Grid[y][x].Alt);
            }

            /// @brief Extract Attributes - Extract a string input into a set of booleans
            /// @param input Set of attributes to unapply (in any order): bo = Bold, it = Italic, un = Underline, re = Reverse, bl = Blink, di = Dim, in = Invisible, st = Standout, pr = Protected, al = Altset
            /// @returns A vector (size of 10) of booleans representing attributes
            std::vector<bool> extractAttributes(std::string input) {
                for (unsigned int i = 0; i < input.length(); i++) {
                    input[i] = tolower(input[i]);
                }

                std::vector<bool> output = {false, false, false, false, false, false, false, false, false, false};

                if (input.find("bo") != std::string::npos) {output[0] = true;}
                if (input.find("it") != std::string::npos) {output[1] = true;}
                if (input.find("un") != std::string::npos) {output[2] = true;}
                if (input.find("re") != std::string::npos) {output[3] = true;}
                if (input.find("bl") != std::string::npos) {output[4] = true;}
                if (input.find("di") != std::string::npos) {output[5] = true;}
                if (input.find("in") != std::string::npos) {output[6] = true;}
                if (input.find("st") != std::string::npos) {output[7] = true;}
                if (input.find("pr") != std::string::npos) {output[8] = true;}
                if (input.find("al") != std::string::npos) {output[9] = true;}

                return output;
            }
            /// @brief Toggle Cell Attributes - Toggle the attributes being used when writing to a cell - boolean based (acts as wattron() or wattroff())
            /// @param bold Bold
            /// @param italic Italic
            /// @param under Underline
            /// @param rev Reverse
            /// @param blink Blink
            /// @param dim Dim
            /// @param invis Invisible
            /// @param stand Standout
            /// @param prot Protected
            /// @param alt Alt Char Set
            void toggleAttributes(bool bold, bool italic, bool under, bool rev, bool blink, bool dim, bool invis, bool stand, bool prot, bool alt) {
                if (bold) {wattron(Win, A_BOLD);}
                else {wattroff(Win, A_BOLD);}
                if (italic) {wattron(Win, A_ITALIC);}
                else {wattroff(Win, A_ITALIC);}
                if (under) {wattron(Win, A_UNDERLINE);}
                else {wattroff(Win, A_UNDERLINE);}
                if (rev) {wattron(Win, A_REVERSE);}
                else {wattroff(Win, A_REVERSE);}
                if (blink) {wattron(Win, A_BLINK);}
                else {wattroff(Win, A_BLINK);}
                if (dim) {wattron(Win, A_DIM);}
                else {wattroff(Win, A_DIM);}
                if (invis) {wattron(Win, A_INVIS);}
                else {wattroff(Win, A_INVIS);}
                if (stand) {wattron(Win, A_STANDOUT);}
                else {wattroff(Win, A_STANDOUT);}
                if (prot) {wattron(Win, A_PROTECT);}
                else {wattroff(Win, A_PROTECT);}
                if (alt) {wattron(Win, A_ALTCHARSET);}
                else {wattroff(Win, A_ALTCHARSET);}
            }
            /// @brief Toggle Cell Attributes - Toggle the attributes being used when writing to a cell - string based (acts as wattron() or wattroff())
            /// @param input A string that should consist of two-letter keys for which attributes to toggle
            void toggleAttributes(std::string input) {
                std::vector<bool> att = extractAttributes(input);
                toggleAttributes(att[0], att[1], att[2], att[3], att[4], att[5], att[6], att[7], att[8], att[9]);
            }

            //
            // TIMING
            //

            /// @brief Wait - Pause everything for an amount of milliseconds - base
            /// @param millis Milliseconds to wait for
            /// @param skippable Whether the wait can be skipped (with a keyboard input) or not
            /// @returns True if the wait was skipped, false if the wait wasn't
            bool wait(unsigned long millis) {
                for (unsigned long i = 0; i < millis; i++) {
                    if (CanSkip) {
                        if (gchar(false, true, false) != ERR) {return true;}
                    }
                    napms(1);
                }

                return false;
            }

            //
            // COMMON CHECKS & CODE SHORTCUTS
            //

            /// @brief Check Coordinate - Check if the coordinates are within the window or not
            /// @param pos Pair consisting of a y-position (row) and an x-position (col) of the coordinate to check
            /// @returns True if the coordinates are in the window, false if the coordinates aren't
            bool checkCoord(std::pair<unsigned short, unsigned short> pos, bool usePadding) {
                if (usePadding) {return !(pos.first < PadUp || pos.first >= DimY - PadDown || pos.second < PadLeft || pos.second >= DimX - PadRight);}
                return !(pos.first < 0 || pos.first >= DimY || pos.second < 0 || pos.second >= DimX);
            }
            /// @brief Check Coordinate - Check if the coordinates are within the window or not
            /// @param y y-position (row) of the coordinate to check
            /// @param x x-position (col) of the coordinate to check
            /// @returns True if the coordinates are in the window, false if the coordinates aren't
            bool checkCoord(unsigned short y, unsigned short x, bool usePadding) {return checkCoord({y, x}, usePadding);}

            //
            // LINE DRAWING HELPERS
            //

            /// @brief Get Connect Style - Get a number corresponding to the type of line that would be needed to connect to the specified character in the specified direction
            /// @param y y-position (row) of the character to get a connection style from
            /// @param x x-position (col) of the character to get a connection style from
            /// @param dir The direction to check: 0 = Up/North, 1 = Down/South, 2 = Left/West, 3 = Right/East
            /// @returns A number 0-3: 0 = Nothing, 1 = Light, 2 = Heavy, 3 = Double
            const unsigned char getConnectStyle(unsigned short y, unsigned short x, unsigned char dir, bool padding) {
                if (!checkCoord(y, x, false) || dir < 0 || dir > 3) {return 0;}

                // If the targeted cell can't be merged with (which is the case for most cells)
                // then there must not be something to connect with (even if the targeted cell is a box drawing element)
                if (!(padding ? smerge(y + PadUp, x + PadLeft) : smerge(y, x))) {return 0;}
                
                wchar_t character = padding ? schar(y + PadUp, x + PadLeft) : schar(y, x);
                for (unsigned char i = 0; i < Keys.LineDrawing.size(); i++) {
                    if (character == Keys.LineDrawing[i].first) {return Keys.LineDrawing[i].second[dir];}
                }

                return 0;
            }

            /// @brief Get Piece - Convert a set of directions (connection styles) and a line style into the corresponding unicode character
            /// @param dir A set of four directions (up, down, left, right)
            /// @param style A pair consisting of a number corresponding to the main style (0 = Light Hard, 1 = Light Soft, 2 = Heavy Both, 3 = Heavy Horizontal, 4 = Heavy Vertical, 5 = Doubled Both, 6 = Doubled Horizontal, 7 = Doubled Vertical) and a number corresponding to the dash style (0 = None, 1 = Double, 2 = Triple, 3 = Quardruple, 4 = Single 1, 5 = Single 2)
            /// @returns A unicode character (box drawing) if a match to the directions is found, a space if a match isn't found
            const wchar_t getPiece(std::vector<unsigned char> dir, std::pair<unsigned char, unsigned char> style) {
                for (unsigned char i = 0; i < Keys.LineDrawing.size(); i++) {
                    if (dir == Keys.LineDrawing[i].second) {
                        // Vertical and horizontal lines (except for the doubled styles)
                        if (i < 24) {
                            if (style.first == HEAVY_HORIZONTAL) {return Keys.LineDrawing[i == 18 ? 12 + style.second : i + style.second].first;}
                            else if (style.first == HEAVY_VERTICAL) {return Keys.LineDrawing[i == 6 ? style.second : i + style.second].first;}
                            
                            return Keys.LineDrawing[i + style.second].first;
                        }
                        // Vertical and horizontal lines for the doubled styles (they don't have dashed variants so are on their own)
                        else if (i < 26) {
                            if (style.first == DOUBLED_HORIZTONAL) {return Keys.LineDrawing[i == 24 ? i : style.second].first;}
                            else if (style.first == DOUBLED_VERTICAL) {return Keys.LineDrawing[i == 24 ? 12 + style.second : i].first;}

                            return Keys.LineDrawing[i].first;
                        }
                        // Corner pieces for the LIGHT_HARD and LIGHT_SOFT styles
                        else if (i < 33) {
                            return Keys.LineDrawing[i + (style.first == LIGHT_SOFT ? 1 : 0)].first;
                        }

                        return Keys.LineDrawing[i].first;
                    }
                }

                return L' ';
            }

        public:
            /// @brief Create an ncursespp Window
            /// @param y Y position (row) for the top-left corner of the Window
            /// @param x X position (col) for the top-left corner of the Window
            /// @param dimy Height (rows) of the Window
            /// @param dimx Length (cols) of the Window
            Window(unsigned short y = 0, unsigned short x = 0, unsigned short dimy = LINES, unsigned short dimx = COLS) {
                // Prevent the window from being made out of bounds (and automatically resize ones that may)
                y = (y < 0 || y > LINES) ? 0 : y;
                x = (x < 0 || x > COLS) ? 0 : x;
                dimy = (dimy < 1 || LINES - dimy - y < 0) ? LINES - y : dimy;
                dimx = (dimx < 1 || COLS - dimx - x < 0) ? COLS - x : dimx;

                PosY = y;
                PosX = x;
                DimY = dimy;
                DimX = dimx;
                Win = newwin(DimY, DimX, PosY, PosX);

                for (unsigned short i = 0; i < DimY; i++) {
                    Grid.emplace_back();
                    for (unsigned short j = 0; j < DimX; j++) {
                        Grid[i].emplace_back();
                    }
                }
            }
            /// @brief Create an auto-centered Window relative to the screen
            /// @param dimy Height (rows) of the Window
            /// @param dimx Length (cols) of the Window
            Window(unsigned short dimy, unsigned short dimx) {
                unsigned short y = LINES / 2 - dimy / 2;
                unsigned short x = COLS / 2 - dimx / 2;
                Window(y, x, dimy, dimx);
            }
            /// @brief Create an auto-centered Window relative to another window
            /// @param w Another window to center the new one to
            /// @param dimy Height (rows) of the Window
            /// @param dimx Length (cols) of the Window
            Window(Window &win, unsigned short dimy, unsigned short dimx) {
                dimy = (win.gdimy() < dimy) ? win.gdimy() : dimy;
                dimx = (win.gdimx() < dimx) ? win.gdimx() : dimx;
                unsigned short y = win.gposy() + win.gdimy() / 2 - dimy / 2;
                unsigned short x = win.gposx() + win.gdimx() / 2 - dimx / 2;
                Window(y, x, dimy, dimx);
            }

            void move(unsigned short y, unsigned short x) {
                PosY = (y < 0) ? 0 : y;
                PosY = (y + DimY > LINES) ? LINES - DimY : y;
                PosX = (x < 0) ? 0 : x;
                PosX = (x + DimX > COLS) ? COLS - DimX : x;

                // Remake the ncurses window
                wclear(Win);
                wrefresh(Win);
                delwin(Win);
                Win = newwin(DimY, DimX, PosY, PosX);
            }

            void resize(unsigned short dimy, unsigned short dimx) {
                dimy = (dimy < 1) ? 1 : dimy;
                dimy = (dimy + PosY > LINES) ? LINES - PosY : dimy;
                dimx = (dimx < 1) ? 1 : dimx;
                dimx = (dimx + PosX > COLS) ? COLS - PosX : dimx;

                std::vector<std::vector<npp::Cell>> buff = Grid;
                Grid.clear();

                for (unsigned short i = 0; i < ((dimy > DimY) ? dimy : DimY); i++) {
                    Grid.emplace_back();
                    for (unsigned short j = 0; j < ((dimx > DimX) ? dimx : DimX); j++) {
                        Grid[i].emplace_back();
                    }
                }
                for (unsigned short i = 0; i < ((dimy > DimY) ? DimY : dimy); i++) {
                    for (unsigned short j = 0; j < ((dimx > DimX) ? DimX : dimx); j++) {
                        Grid[i][j] = buff[i][j];
                    }
                }

                DimY = dimy;
                DimX = dimx;

                // Remake the ncurses window
                wclear(Win);
                wrefresh(Win);
                delwin(Win);
                Win = newwin(DimY, DimX, PosY, PosX);
            }

            //
            // MISC
            //

            /// @brief Return Wait - Pause everything for an amount of milliseconds - char return
            /// @param millis Milliseconds to wait for
            /// @returns An integer that can be used in the same way that it is in the ncurses getch() function
            int rwait(unsigned long millis) {
                int input;

                for (unsigned long i = 0; i < millis; i++) {
                    if (CanSkip) {
                        input = gchar(false, true, false);

                        if (input != ERR) {return input;}
                    }
                    napms(1);
                }

                return -1;
            }
            /// @brief Target Wait - Pause everything for an amount of milliseconds - target
            /// @param millis Milliseconds to wait for
            /// @param targets A list of characters that would skip the wait if detected
            /// @returns True if the wait was skipped, false if the wait wasn't
            bool twait(unsigned long millis, std::vector<char> targets) {
                int input;

                for (unsigned long i = 0; i < millis; i++) {
                    if (CanSkip) {
                        input = gchar(false, true, false);

                        if (input != ERR) {
                            for (unsigned short j = 0; j < targets.size(); j++) {
                                if (input == targets[j]) {return true;}
                            }
                        }
                    }
                    napms(1);
                }

                return false;
            }

            /// @brief Update Skippability - Change whether the window will allow the user to skip wait() functions
            /// @param skippable Whether the window will allow the user to skip or not
            void uskip(bool skippable) {CanSkip = skippable;}

            /// @brief Update Padding (Top) - Change how many padding spaces will be on the top of the window
            /// @param padding Amount of characters to pad
            void upup(unsigned short padding = 0) {PadUp = padding < 0 ? PadUp : padding;}
            /// @brief Update Padding (Bottom) - Change how many padding spaces will be on the bottom of the window
            /// @param padding Amount of characters to pad
            void updown(unsigned short padding = 0) {PadDown = padding < 0 ? PadDown : padding;}
            /// @brief Update Padding (Left) - Change how many padding spaces will be on the left of the window
            /// @param padding Amount of characters to pad
            void upleft(unsigned short padding = 0) {PadLeft = padding < 0 ? PadLeft : padding;}
            /// @brief Update Padding (Right) - Change how many padding spaces will be on the right of the window
            /// @param padding Amount of characters to pad
            void upright(unsigned short padding = 0) {PadRight = padding < 0 ? PadRight : padding;}
            /// @brief Update Padding (All) - Change how many padding spaces will be on all sides of the window
            /// @param padding 
            void upall(unsigned short padding = 0) {if (padding >= 0) {PadUp = PadDown = PadLeft = PadRight = padding;}}

            //
            // GETTING WINDOW/CELL ATTRIBUTES
            //

            /// @brief Get Y-Dimension - Get the y-dimension (rows) of a window
            /// @returns The y-dimension (rows) of a window
            const unsigned short gdimy(bool padding = Defaults.UsePadding) {return DimY - (padding ? PadUp + PadDown : 0);}
            /// @brief Get X-Dimension - Get the x-dimension (cols) of a window
            /// @returns The x-dimension (cols) of a window
            const unsigned short gdimx(bool padding = Defaults.UsePadding) {return DimX - (padding ? PadLeft + PadRight : 0);}
            /// @brief Get Y-Position - Get the y-position (row) of the top-left corner of a window
            /// @returns The y-dimension (row) of the top-left corner of a window
            const unsigned short gposy() {return PosY;}
            /// @brief Get X-Position - Get the x-position (col) of the top-left corner of a window
            /// @returns The x-dimension (col) of the top-left corner of a window
            const unsigned short gposx() {return PosX;}
            /// @brief Get Padding (Top) - Get the amount of characters padding the top of the screen
            /// @returns The amount of characters padding the top of the screen
            const unsigned short gpadt() {return PadUp;}
            /// @brief Get Padding (Bottom) - Get the amount of characters padding the bottom of the screen
            /// @returns The amount of characters padding the bottom of the screen
            const unsigned short gpadb() {return PadDown;}
            /// @brief Get Padding (Left) - Get the amount of characters padding the left of the screen
            /// @returns The amount of characters padding the left of the screen
            const unsigned short gpadl() {return PadLeft;}
            /// @brief Get Padding (Right) - Get the amount of characters padding the right of the screen
            /// @returns The amount of characters padding the right of the screen
            const unsigned short gpadr() {return PadRight;}

            /// @brief Scan Character - Get the character of a cell from a window
            /// @param y y-position (row) of the cell to scan
            /// @param x x-position (col) of the cell to scan
            /// @returns A string (length of 1) representing a cell's character or a space if the requested cell doesn't exist
            const wchar_t schar(unsigned short y, unsigned short x) {return !checkCoord(y, x, false) ? L' ' : Grid[y][x].Char;}
            
            /// @brief Scan Color - Get the color pair of a cell from a window
            /// @param y y-position (row) of the cell to scan
            /// @param x x-position (col) of the cell to scan
            /// @returns The number of the color pair a cell uses or the default pair for the window if the requested cell doesn't exist
            const unsigned char scolor(unsigned short y, unsigned short x) {return !checkCoord(y, x, false) ? 1 : Grid[y][x].Color;}
            
            /// @brief Scan Bold - Get whether a cell from a window is bolded or not
            /// @param y y-position (row) of the cell to scan
            /// @param x x-position (col) of the cell to scan
            /// @returns A true or false for if the cell is bolded or the default state for the window if the requested cell doesn't exist
            const bool sbold(unsigned short y, unsigned short x) {return !checkCoord(y, x, false) ? false : Grid[y][x].Bold;}
            
            /// @brief Scan Italic - Get whether a cell from a window is italicized or not
            /// @param y y-position (row) of the cell to scan
            /// @param x x-position (col) of the cell to scan
            /// @returns A true or false for if the cell is italicized or the default state for the window if the requested cell doesn't exist
            const bool sitalic(unsigned short y, unsigned short x) {return !checkCoord(y, x, false) ? false : Grid[y][x].Italic;}
            
            /// @brief Scan Underline - Get whether a cell from a window is underlined or not
            /// @param y y-position (row) of the cell to scan
            /// @param x x-position (col) of the cell to scan
            /// @returns A true or false for if the cell is underlined or the default state for the window if the requested cell doesn't exist
            const bool sunder(unsigned short y, unsigned short x) {return !checkCoord(y, x, false) ? false : Grid[y][x].Under;}
            
            /// @brief Scan Reverse - Get whether a cell from a window is reversed or not
            /// @param y y-position (row) of the cell to scan
            /// @param x x-position (col) of the cell to scan
            /// @returns A true or false for if the cell is reversed or the default state for the window if the requested cell doesn't exist
            const bool srev(unsigned short y, unsigned short x) {return !checkCoord(y, x, false) ? false : Grid[y][x].Rev;}
            
            /// @brief Scan Blink - Get whether a cell from a window is blinking or not
            /// @param y y-position (row) of the cell to scan
            /// @param x x-position (col) of the cell to scan
            /// @returns A true or false for if the cell is blinking or the default state for the window if the requested cell doesn't exist
            const bool sblink(unsigned short y, unsigned short x) {return !checkCoord(y, x, false) ? false : Grid[y][x].Blink;}
            
            /// @brief Scan Dim - Get whether a cell from a window is dim or not
            /// @param y y-position (row) of the cell to scan
            /// @param x x-position (col) of the cell to scan
            /// @returns A true or false for if the cell is dim or the default state for the window if the requested cell doesn't exist
            const bool sdim(unsigned short y, unsigned short x) {return !checkCoord(y, x, false) ? false : Grid[y][x].Dim;}
            
            /// @brief Scan Invisible - Get whether a cell from a window is invisible or not
            /// @param y y-position (row) of the cell to scan
            /// @param x x-position (col) of the cell to scan
            /// @returns A true or false for if the cell is invisible or the default state for the window if the requested cell doesn't exist
            const bool sinvis(unsigned short y, unsigned short x) {return !checkCoord(y, x, false) ? false : Grid[y][x].Invis;}
            
            /// @brief Scan Standout - Get whether a cell from a window is standing out or not
            /// @param y y-position (row) of the cell to scan
            /// @param x x-position (col) of the cell to scan
            /// @returns A true or false for if the cell is standing out or the default state for the window if the requested cell doesn't exist
            const bool sstand(unsigned short y, unsigned short x) {return !checkCoord(y, x, false) ? false : Grid[y][x].Stand;}
            
            /// @brief Scan Protected - Get whether a cell from a window is protected or not
            /// @param y y-position (row) of the cell to scan
            /// @param x x-position (col) of the cell to scan
            /// @returns A true or false for if the cell is protected or the default state for the window if the requested cell doesn't exist
            const bool sprot(unsigned short y, unsigned short x) {return !checkCoord(y, x, false) ? false : Grid[y][x].Prot;}
            
            /// @brief Scan Alt Char Set - Get whether a cell from a window is using the alternate character set or not
            /// @param y y-position (row) of the cell to scan
            /// @param x x-position (col) of the cell to scan
            /// @returns A true or false for if the cell is using the alternate character set or the default state for the window if the requested cell doesn't exist
            const bool salt(unsigned short y, unsigned short x) {return  !checkCoord(y, x, false) ? false : Grid[y][x].Alt;}

            /// @brief Scan Mergeability - Get whether a cell can merge with other cells (should only apply to cells that contain box drawing characters)
            /// @param y y-position (row) of the cell to scan
            /// @param x x-position (col) of the cell to scan
            /// @returns A true or false for if the cell can merge with other cells
            const bool smerge(unsigned short y, unsigned short x) {return !checkCoord(y, x, false) ? false : Grid[y][x].CanMerge;}

            const npp::Cell scell(unsigned short y, unsigned short x) {return !checkCoord(y, x, false) ? npp::Cell() : Grid[y][x];}

            //
            // WRITING TO WINDOW
            //

            /// @brief Clear the window without removing cell data
            void clear() {wclear(Win);}
            /// @brief Clear a portion of the window without removing cell data
            /// @param y y-position (row) to start clearing the window at
            /// @param x x-position (col) to start clearing the window at
            /// @param dimy Amount of rows to clear
            /// @param dimx Amount of columns to clear
            /// @param usePadding Whether to offset coordinates and disallow the editing of certain cells or not
            void clear(unsigned short y, unsigned short x, unsigned short dimy, unsigned short dimx, bool usePadding = Defaults.UsePadding) {
                if (!checkCoord(y, x, usePadding) || !checkCoord(y + dimy - 1, x + dimx - 1, usePadding)) {return;}

                for (unsigned short i = 0; i < dimy; i++) {
                    for (unsigned short j = 0; j < dimx; j++) {
                        mvwaddch(Win, y + i, x + j, ' ');
                    }
                }
            }
            /// @brief Clear the window and remove cell data
            void reset() {
                for (unsigned short i = 0; i < DimY; i++) {
                    for (unsigned short j = 0; j < DimX; j++) {
                        Grid[i][j] = Cell();
                    }
                }
                clear();
            }
            /// @brief Clear a portion of the window and remove cell data
            /// @param y y-position (row) to start reseting the window at
            /// @param x x-position (col) to start reseting the window at
            /// @param dimy Amount of rows to reset
            /// @param dimx Amount of columns to reset
            /// @param usePadding Whether to offset coordinates and disallow the editing of certain cells or not
            void reset(unsigned short y, unsigned short x, unsigned short dimy, unsigned short dimx, bool usePadding = Defaults.UsePadding) {
                if (!checkCoord(y, x, usePadding) || !checkCoord(y + dimy - 1, x + dimx - 1, usePadding)) {return;}

                for (unsigned short i = 0; i < dimy; i++) {
                    for (unsigned short j = 0; j < dimx; j++) {
                        Grid[y + i][x + j] = Cell();
                    }
                }
                clear(y, x, dimy, dimx, usePadding);
            }

            /// @brief Write Character, Return Position - Write a single character to the window - pair pos, pair return
            /// @param pos Pair consisting of a y-position (row) and an x-position (col) for the character to be written at
            /// @param input Wide character input to be written (unicode friendly)
            /// @param color Color pair to use when writing each cell
            /// @param att Set of attributes to apply (in any order): bo = Bold, it = Italic, un = Underline, re = Reverse, bl = Blink, di = Dim, in = Invisible, st = Standout, pr = Protected, al = Altset
            /// @param offset Pair consisiting of a y-offset (rows) and an x-offset (cols) that changes the return value
            /// @param usePadding Whether to disallow writing to cells that are in a padded region or not
            /// @returns A pair consisting of a y-position (row) and an x-position (col) based off of the character's position and the inputted offset
            std::pair<unsigned short, unsigned short> wcharp(std::pair<unsigned short, unsigned short> pos, wchar_t input, unsigned char color = Defaults.Color, std::string att = Defaults.Attributes, std::pair<unsigned short, unsigned short> offset = Defaults.Offset, bool padding = Defaults.UsePadding) {
                if (padding) {
                    pos.first += PadUp;
                    pos.second += PadLeft;
                    if (!checkCoord(pos, true)) {return {pos.first - PadUp, pos.second - PadLeft};}
                } else {
                    if (!checkCoord(pos, false)) {return {pos.first, pos.second};}
                }

                std::vector<bool> attributes = extractAttributes(att);

                Grid[pos.first][pos.second].Char = input;
                Grid[pos.first][pos.second].Color = color;
                
                Grid[pos.first][pos.second].Bold = attributes[0];
                Grid[pos.first][pos.second].Italic = attributes[1];
                Grid[pos.first][pos.second].Under = attributes[2];
                Grid[pos.first][pos.second].Rev = attributes[3];
                Grid[pos.first][pos.second].Blink = attributes[4];
                Grid[pos.first][pos.second].Dim = attributes[5];
                Grid[pos.first][pos.second].Invis = attributes[6];
                Grid[pos.first][pos.second].Stand = attributes[7];
                Grid[pos.first][pos.second].Prot = attributes[8];
                Grid[pos.first][pos.second].Alt = attributes[9];

                Grid[pos.first][pos.second].CanMerge = false;

                if (padding) {return {pos.first + offset.first - PadUp, pos.second + offset.second - PadLeft};}
                return {pos.first + offset.first, pos.second + offset.second};
            }
            /// @brief Write Character, Return Position - Write a single character (wchar_t) to the window - short short pos, pair return
            /// @param y y-position (row) of the character
            /// @param x x-position (col) of the character
            /// @param input Wide character input to be written (unicode friendly)
            /// @param color Color pair to use when writing each cell
            /// @param att Set of attributes to apply (in any order): bo = Bold, it = Italic, un = Underline, re = Reverse, bl = Blink, di = Dim, in = Invisible, st = Standout, pr = Protected, al = Altset
            /// @param offset Pair consisiting of a y-offset (rows) and an x-offset (cols) that changes the return value
            /// @param usePadding Whether to disallow writing to cells that are in a padded region or not
            /// @returns A pair consisting of a y-position (row) and an x-position (col) based off of the character's position and the inputted offset
            std::pair<unsigned short, unsigned short> wcharp(unsigned short y, unsigned short x, wchar_t input, unsigned char color = Defaults.Color, std::string att = Defaults.Attributes, std::pair<unsigned short, unsigned short> offset = Defaults.Offset, bool padding = Defaults.UsePadding) {return wcharp({y, x}, input, color, att, offset, padding);}
            /// @brief Write Character, Return Y-Position - Write a single character to the window - pair pos, short return
            /// @param pos Pair consisting of a y-position (row) and an x-position (col) for the character to be written at
            /// @param input Wide character input to be written (unicode friendly)
            /// @param color Color pair to use when writing each cell
            /// @param att Set of attributes to apply (in any order): bo = Bold, it = Italic, un = Underline, re = Reverse, bl = Blink, di = Dim, in = Invisible, st = Standout, pr = Protected, al = Altset
            /// @param offset Pair consisiting of a y-offset (rows) and an x-offset (cols) that changes the return value
            /// @param usePadding Whether to disallow writing to cells that are in a padded region or not
            /// @returns A y-position (row) based off of the character's position and the inputted offset
            unsigned short wchary(std::pair<unsigned short, unsigned short> pos, wchar_t input, unsigned char color = Defaults.Color, std::string att = Defaults.Attributes, std::pair<unsigned short, unsigned short> offset = Defaults.Offset, bool padding = Defaults.UsePadding) {return wcharp(pos, input, color, att, offset, padding).first;}
            /// @brief Write Character, Return Y-Position - Write a single character to the window - short short pos, short return
            /// @param y y-position (row) of the character
            /// @param x x-position (col) of the character
            /// @param input Wide character input to be written (unicode friendly)
            /// @param color Color pair to use when writing each cell
            /// @param att Set of attributes to apply (in any order): bo = Bold, it = Italic, un = Underline, re = Reverse, bl = Blink, di = Dim, in = Invisible, st = Standout, pr = Protected, al = Altset
            /// @param offset Pair consisiting of a y-offset (rows) and an x-offset (cols) that changes the return value
            /// @param usePadding Whether to disallow writing to cells that are in a padded region or not
            /// @returns A y-position (row) based off of the character's position and the inputted offset
            unsigned short wchary(unsigned short y, unsigned short x, wchar_t input, unsigned char color = Defaults.Color, std::string att = Defaults.Attributes, std::pair<unsigned short, unsigned short> offset = Defaults.Offset, bool padding = Defaults.UsePadding) {return wchary({y, x}, input, color, att, offset, padding);}
            /// @brief Write Character, Return X-Position - Write a single character to the window - pair pos, short return
            /// @param pos Pair consisting of a y-position (row) and an x-position (col) for the character to be written at
            /// @param input Wide character input to be written (unicode friendly)
            /// @param color Color pair to use when writing each cell
            /// @param att Set of attributes to apply (in any order): bo = Bold, it = Italic, un = Underline, re = Reverse, bl = Blink, di = Dim, in = Invisible, st = Standout, pr = Protected, al = Altset
            /// @param offset Pair consisiting of a y-offset (rows) and an x-offset (cols) that changes the return value
            /// @param usePadding Whether to disallow writing to cells that are in a padded region or not
            /// @returns A x-position (col) based off of the character's position and the inputted offset
            unsigned short wcharx(std::pair<unsigned short, unsigned short> pos, wchar_t input, unsigned char color = Defaults.Color, std::string att = Defaults.Attributes, std::pair<unsigned short, unsigned short> offset = Defaults.Offset, bool padding = Defaults.UsePadding) {return wcharp(pos, input, color, att, offset, padding).second;}
            /// @brief Write Character, Return X-Position - Write a single character to the window - short short pos, short return
            /// @param y y-position (row) of the character
            /// @param x x-position (col) of the character
            /// @param input Wide character input to be written (unicode friendly)
            /// @param color Color pair to use when writing each cell
            /// @param att Set of attributes to apply (in any order): bo = Bold, it = Italic, un = Underline, re = Reverse, bl = Blink, di = Dim, in = Invisible, st = Standout, pr = Protected, al = Altset
            /// @param offset Pair consisiting of a y-offset (rows) and an x-offset (cols) that changes the return value
            /// @param usePadding Whether to disallow writing to cells that are in a padded region or not
            /// @returns A x-position (col) based off of the character's position and the inputted offset
            unsigned short wcharx(unsigned short y, unsigned short x, wchar_t input, unsigned char color = Defaults.Color, std::string att = Defaults.Attributes, std::pair<unsigned short, unsigned short> offset = Defaults.Offset, bool padding = Defaults.UsePadding) {return wcharx({y, x}, input, color, att, offset, padding);}
            /// @brief Write Character, Return Nothing - Write a single character to the window - pair pos, no return
            /// @param pos Pair consisting of a y-position (row) and an x-position (col) for the character to be written at
            /// @param input Wide character input to be written (unicode friendly)
            /// @param color Color pair to use when writing each cell
            /// @param att Set of attributes to apply (in any order): bo = Bold, it = Italic, un = Underline, re = Reverse, bl = Blink, di = Dim, in = Invisible, st = Standout, pr = Protected, al = Altset
            /// @param offset Pair consisiting of a y-offset (rows) and an x-offset (cols) that changes the return value
            /// @param usePadding Whether to disallow writing to cells that are in a padded region or not
            void wchar(std::pair<unsigned short, unsigned short> pos, wchar_t input, unsigned char color = Defaults.Color, std::string att = Defaults.Attributes, bool padding = Defaults.UsePadding) {wcharp(pos, input, color, att, {0, 0}, padding);}
            /// @brief Write Character, Return Nothing - Write a single character to the window - short short pos, no return
            /// @param y y-position (row) of the character
            /// @param x x-position (col) of the character
            /// @param input Wide character input to be written (unicode friendly)
            /// @param color Color pair to use when writing each cell
            /// @param att Set of attributes to apply (in any order): bo = Bold, it = Italic, un = Underline, re = Reverse, bl = Blink, di = Dim, in = Invisible, st = Standout, pr = Protected, al = Altset
            /// @param offset Pair consisiting of a y-offset (rows) and an x-offset (cols) that changes the return value
            /// @param usePadding Whether to disallow writing to cells that are in a padded region or not
            void wchar(unsigned short y, unsigned short x, wchar_t input, unsigned char color = Defaults.Color, std::string att = Defaults.Attributes, bool padding = Defaults.UsePadding) {wchar({y, x}, input, color, att, padding);}

            /// @brief Write String, Return Position - Write a string to the window - pair pos, pair return
            /// @param pos Pair consisting of a y-position (row) and an x-position (col) indicating the start of the string
            /// @param input String input to be written (unicode friendly)
            /// @param color Color pair to use when writing each cell
            /// @param att Set of attributes to apply (in any order): bo = Bold, it = Italic, un = Underline, re = Reverse, bl = Blink, di = Dim, in = Invisible, st = Standout, pr = Protected, al = Altset
            /// @param offset Pair consisiting of a y-offset (rows) and an x-offset (cols) that changes the return value
            /// @param wrap Whether to wrap to starting x or not
            /// @param usePadding Whether to disallow writing to cells that are in a padded region or not
            /// @returns A pair consisting of a y-position (row) and an x-position (col) based off of the last character of the string that was written and the inputted offset
            std::pair<unsigned short, unsigned short> wstrp(std::pair<unsigned short, unsigned short> pos, std::wstring input, unsigned char color = Defaults.Color, std::string att = Defaults.Attributes, std::pair<unsigned short, unsigned short> offset = Defaults.Offset, bool wrap = Defaults.Wrap, bool padding = Defaults.UsePadding) {
                if (!checkCoord(pos, false)) {return {pos.first, pos.second};}

                unsigned short yoffset = 0, xoffset = 0;
                for (unsigned short i = 0; i < input.length(); i++) {
                    // Automatically wrap the string back to the starting x-pos
                    if (pos.second + i - xoffset >= DimX - (padding ? PadLeft + PadRight : 0)) {
                        yoffset++;
                        xoffset += DimX - (wrap ? pos.second : 0) - (padding ? PadLeft + PadRight : 0);
                    }

                    // Quit early if the string goes out of bounds
                    if (!checkCoord(pos.first + yoffset, pos.second + i - xoffset, false)) {break;}

                    wchar(pos.first + yoffset, pos.second + i - xoffset, input[i], color, att, padding);
                }

                return {pos.first + yoffset + offset.first, pos.second + input.length() - 1 - xoffset + offset.second};
            }
            /// @brief Write String, Return Position - Write a string to the window - short short pos, pair return
            /// @param y y-position (row) of the start of the string
            /// @param x x-position (col) of the start of the string
            /// @param input String input to be written (unicode friendly)
            /// @param color Color pair to use when writing each cell
            /// @param att Set of attributes to unapply (in any order): bo = Bold, it = Italic, un = Underline, re = Reverse, bl = Blink, di = Dim, in = Invisible, st = Standout, pr = Protected, al = Altset
            /// @param offset Pair consisiting of a y-offset (rows) and an x-offset (cols) that changes the return value
            /// @param wrap Whether to wrap to starting x or not
            /// @param usePadding Whether to disallow writing to cells that are in a padded region or not
            /// @returns A pair consisting of a y-position (row) and an x-position (col) based off of the last character of the string that was written and the inputted offset
            std::pair<unsigned short, unsigned short> wstrp(unsigned short y, unsigned short x, std::wstring input, unsigned char color = Defaults.Color, std::string att = Defaults.Attributes, std::pair<unsigned short, unsigned short> offset = Defaults.Offset, bool wrap = Defaults.Wrap, bool padding = Defaults.UsePadding) {return wstrp({y, x}, input, color, att, offset, wrap, padding);}
            /// @brief Write String, Return Y-Position - Write a string to the window - pair pos, short return
            /// @param pos Pair consisting of a y-position (row) and an x-position (col) indicating the start of the string
            /// @param input String input to be written (unicode friendly)
            /// @param color Color pair to use when writing each cell
            /// @param att Set of attributes to apply (in any order): bo = Bold, it = Italic, un = Underline, re = Reverse, bl = Blink, di = Dim, in = Invisible, st = Standout, pr = Protected, al = Altset
            /// @param offset Pair consisiting of a y-offset (rows) and an x-offset (cols) that changes the return value
            /// @param wrap Whether to wrap to starting x or not
            /// @param usePadding Whether to disallow writing to cells that are in a padded region or not
            /// @returns A y-position (row) based off of the last character of the string that was written and the inputted offset
            unsigned short wstry(std::pair<unsigned short, unsigned short> pos, std::wstring input, unsigned char color = Defaults.Color, std::string att = Defaults.Attributes, std::pair<unsigned short, unsigned short> offset = Defaults.Offset, bool wrap = Defaults.Wrap, bool padding = Defaults.UsePadding) {return wstrp(pos, input, color, att, offset, wrap, padding).first;}
            /// @brief Write String, Return Y-Position - Write a string to the window - short short pos, short return
            /// @param y y-position (row) of the start of the string
            /// @param x x-position (col) of the start of the string
            /// @param input String input to be written (unicode friendly)
            /// @param color Color pair to use when writing each cell
            /// @param att Set of attributes to apply (in any order): bo = Bold, it = Italic, un = Underline, re = Reverse, bl = Blink, di = Dim, in = Invisible, st = Standout, pr = Protected, al = Altset
            /// @param offset Pair consisiting of a y-offset (rows) and an x-offset (cols) that changes the return value
            /// @param wrap Whether to wrap to starting x or not
            /// @param usePadding Whether to disallow writing to cells that are in a padded region or not
            /// @returns A y-position (row) based off of the last character of the string that was written and the inputted offset
            unsigned short wstry(unsigned short y, unsigned short x, std::wstring input, unsigned char color = Defaults.Color, std::string att = Defaults.Attributes, std::pair<unsigned short, unsigned short> offset = Defaults.Offset, bool wrap = Defaults.Wrap, bool padding = Defaults.UsePadding) {return wstry({y, x}, input, color, att, offset, wrap, padding);}
            /// @brief Write String, Return X-Position - Write a string to the window - pair pos, short return
            /// @param pos Pair consisting of a y-position (row) and an x-position (col) indicating the start of the string
            /// @param input String input to be written (unicode friendly)
            /// @param color Color pair to use when writing each cell
            /// @param att Set of attributes to apply (in any order): bo = Bold, it = Italic, un = Underline, re = Reverse, bl = Blink, di = Dim, in = Invisible, st = Standout, pr = Protected, al = Altset
            /// @param offset Pair consisiting of a y-offset (rows) and an x-offset (cols) that changes the return value
            /// @param wrap Whether to wrap to starting x or not
            /// @param usePadding Whether to disallow writing to cells that are in a padded region or not
            /// @returns A x-position (col) based off of the last character of the string that was written and the inputted offset
            unsigned short wstrx(std::pair<unsigned short, unsigned short> pos, std::wstring input, unsigned char color = Defaults.Color, std::string att = Defaults.Attributes, std::pair<unsigned short, unsigned short> offset = Defaults.Offset, bool wrap = Defaults.Wrap, bool padding = Defaults.UsePadding) {return wstrp(pos, input, color, att, offset, wrap, padding).second;}
            /// @brief Write String, Return X-Position - Write a string to the window - short short pos, short return
            /// @param y y-position (row) of the start of the string
            /// @param x x-position (col) of the start of the string
            /// @param input String input to be written (unicode friendly)
            /// @param color Color pair to use when writing each cell
            /// @param att Set of attributes to apply (in any order): bo = Bold, it = Italic, un = Underline, re = Reverse, bl = Blink, di = Dim, in = Invisible, st = Standout, pr = Protected, al = Altset
            /// @param offset Pair consisiting of a y-offset (rows) and an x-offset (cols) that changes the return value
            /// @param wrap Whether to wrap to starting x or not
            /// @param usePadding Whether to disallow writing to cells that are in a padded region or not
            /// @returns A x-position (col) based off of the last character of the string that was written and the inputted offset
            unsigned short wstrx(unsigned short y, unsigned short x, std::wstring input, unsigned char color = Defaults.Color, std::string att = Defaults.Attributes, std::pair<unsigned short, unsigned short> offset = Defaults.Offset, bool wrap = Defaults.Wrap, bool padding = Defaults.UsePadding) {return wstrx({y, x}, input, color, att, offset, wrap, padding);}
            /// @brief Write String, Return Nothing - Write a string to the window - pair pos, no return
            /// @param pos Pair consisting of a y-position (row) and an x-position (col) indicating the start of the string
            /// @param input String input to be written (unicode friendly)
            /// @param color Color pair to use when writing each cell
            /// @param att Set of attributes to apply (in any order): bo = Bold, it = Italic, un = Underline, re = Reverse, bl = Blink, di = Dim, in = Invisible, st = Standout, pr = Protected, al = Altset
            /// @param wrap Whether to wrap to starting x or not
            /// @param usePadding Whether to disallow writing to cells that are in a padded region or not
            void wstr(std::pair<unsigned short, unsigned short> pos, std::wstring input, unsigned char color = Defaults.Color, std::string att = Defaults.Attributes, bool wrap = Defaults.Wrap, bool padding = Defaults.UsePadding) {wstrp(pos, input, color, att, {0, 0}, wrap, padding);}
            /// @brief Write String, Return Nothing - Write a string to the window - short short pos, no return
            /// @param y y-position (row) of the start of the string
            /// @param x x-position (col) of the start of the string
            /// @param input String input to be written (unicode friendly)
            /// @param color Color pair to use when writing each cell
            /// @param att Set of attributes to apply (in any order): bo = Bold, it = Italic, un = Underline, re = Reverse, bl = Blink, di = Dim, in = Invisible, st = Standout, pr = Protected, al = Altset
            /// @param wrap Whether to wrap to starting x or not
            /// @param usePadding Whether to disallow writing to cells that are in a padded region or not
            void wstr(unsigned short y, unsigned short x, std::wstring input, unsigned char color = Defaults.Color, std::string att = Defaults.Attributes, bool wrap = Defaults.Wrap, bool padding = Defaults.UsePadding) {wstr({y, x}, input, color, att, wrap, padding);}

            /// @brief Write Integer, Return Position - Write an integer to the window - pair pos, pair return
            /// @param pos Pair consisting of a y-position (row) and an x-position (col) indicating the start of the integer
            /// @param input Integer input to be written
            /// @param minWidth Minimum width for the integer to be (leading zeros will be added to meet the width)
            /// @param color Color pair to use when writing each cell
            /// @param att Set of attributes to apply (in any order): bo = Bold, it = Italic, un = Underline, re = Reverse, bl = Blink, di = Dim, in = Invisible, st = Standout, pr = Protected, al = Altset
            /// @param offset Pair consisiting of a y-offset (rows) and an x-offset (cols) that changes the return value
            /// @param wrap Whether to wrap to starting x or not
            /// @param usePadding Whether to disallow writing to cells that are in a padded region or not
            /// @returns A pair consisting of a y-position (row) and x-position (col) of the last character of the integer that was written
            std::pair<unsigned short, unsigned short> wintp(std::pair<unsigned short, unsigned short> pos, int input, unsigned short minWidth = Defaults.MinWidth, unsigned char color = Defaults.Color, std::string att = Defaults.Attributes, std::pair<unsigned short, unsigned short> offset = Defaults.Offset, bool wrap = Defaults.Wrap, bool padding = Defaults.UsePadding) {
                if (!checkCoord(pos, false)) {return pos;}

                // Convert number to a string that can be written to the window
                std::wstring num = std::to_wstring(input);
                unsigned short length = num.length();

                // Add leading zeros; negative numbers are accounted for by keeping the negative sign in front
                for (int i = 0; i < minWidth - length; i++) {
                    num.insert(((input < 0) ? num.begin() + 1 : num.begin()), 1, '0');
                }

                return wstrp({pos.first, pos.second}, num, color, att, offset, wrap, padding);
            }
            /// @brief Write Integer, Return Position - Write an integer to the window - short short pos, pair return
            /// @param y y-position (row) of the start of the integer
            /// @param x x-position (col) of the start of the integer
            /// @param input Integer input to be written
            /// @param minWidth Minimum width for the integer to be (leading zeros will be added to meet the width)
            /// @param color Color pair to use when writing each cell
            /// @param att Set of attributes to apply (in any order): bo = Bold, it = Italic, un = Underline, re = Reverse, bl = Blink, di = Dim, in = Invisible, st = Standout, pr = Protected, al = Altset
            /// @param offset Pair consisiting of a y-offset (rows) and an x-offset (cols) that changes the return value
            /// @param wrap Whether to wrap to starting x or not
            /// @param usePadding Whether to disallow writing to cells that are in a padded region or not
            /// @returns A pair consisting of the y-position (row) and x-position (col) of the last character of the integer that was written
            std::pair<unsigned short, unsigned short> wintp(unsigned short y, unsigned short x, int input, unsigned short minWidth = Defaults.MinWidth, unsigned char color = Defaults.Color, std::string att = Defaults.Attributes, std::pair<unsigned short, unsigned short> offset = Defaults.Offset, bool wrap = Defaults.Wrap, bool padding = Defaults.UsePadding) {return wintp({y, x}, input, minWidth, color, att, offset, wrap, padding);}
            /// @brief Write Integer, Return Y-Position - Write an integer to the window - pair pos, short return
            /// @param pos Pair consisting of a y-position (row) and an x-position (col) indicating the start of the integer
            /// @param input Integer input to be written
            /// @param minWidth Minimum width for the integer to be (leading zeros will be added to meet the width)
            /// @param color Color pair to use when writing each cell
            /// @param att Set of attributes to apply (in any order): bo = Bold, it = Italic, un = Underline, re = Reverse, bl = Blink, di = Dim, in = Invisible, st = Standout, pr = Protected, al = Altset
            /// @param offset Pair consisiting of a y-offset (rows) and an x-offset (cols) that changes the return value
            /// @param wrap Whether to wrap to starting x or not
            /// @param usePadding Whether to disallow writing to cells that are in a padded region or not
            /// @returns A y-position (row) based off of the last character of the integer that was written and the inputted offset
            unsigned short winty(std::pair<unsigned short, unsigned short> pos, int input, unsigned short minWidth = Defaults.MinWidth, unsigned char color = Defaults.Color, std::string att = Defaults.Attributes, std::pair<unsigned short, unsigned short> offset = Defaults.Offset, bool wrap = Defaults.Wrap, bool padding = Defaults.UsePadding) {return wintp(pos, input, minWidth, color, att, offset, wrap, padding).first;}
            /// @brief Write Integer, Return Y-Position - Write an integer to the window - short short pos, short return
            /// @param y y-position (row) of the start of the integer
            /// @param x x-position (col) of the start of the integer
            /// @param input Integer input to be written
            /// @param minWidth Minimum width for the integer to be (leading zeros will be added to meet the width)
            /// @param color Color pair to use when writing each cell
            /// @param att Set of attributes to apply (in any order): bo = Bold, it = Italic, un = Underline, re = Reverse, bl = Blink, di = Dim, in = Invisible, st = Standout, pr = Protected, al = Altset
            /// @param offset Pair consisiting of a y-offset (rows) and an x-offset (cols) that changes the return value
            /// @param wrap Whether to wrap to starting x or not
            /// @param usePadding Whether to disallow writing to cells that are in a padded region or not
            /// @returns A y-position (row) based off of the last character of the integer that was written and the inputted offset
            unsigned short winty(unsigned short y, unsigned short x, int input, unsigned short minWidth = Defaults.MinWidth, unsigned char color = Defaults.Color, std::string att = Defaults.Attributes, std::pair<unsigned short, unsigned short> offset = Defaults.Offset, bool wrap = Defaults.Wrap, bool padding = Defaults.UsePadding) {return winty({y, x}, input, minWidth, color, att, offset, wrap, padding);}
            /// @brief Write Integer, Return X-Position - Write an integer to the window - pair pos, short return
            /// @param pos Pair consisting of a y-position (row) and an x-position (col) indicating the start of the integer
            /// @param input Integer input to be written
            /// @param minWidth Minimum width for the integer to be (leading zeros will be added to meet the width)
            /// @param color Color pair to use when writing each cell
            /// @param att Set of attributes to apply (in any order): bo = Bold, it = Italic, un = Underline, re = Reverse, bl = Blink, di = Dim, in = Invisible, st = Standout, pr = Protected, al = Altset
            /// @param offset Pair consisiting of a y-offset (rows) and an x-offset (cols) that changes the return value
            /// @param wrap Whether to wrap to starting x or not
            /// @param usePadding Whether to disallow writing to cells that are in a padded region or not
            /// @returns A x-position (col) based off of the last character of the integer that was written and the inputted offset
            unsigned short wintx(std::pair<unsigned short, unsigned short> pos, int input, unsigned short minWidth = Defaults.MinWidth, unsigned char color = Defaults.Color, std::string att = Defaults.Attributes, std::pair<unsigned short, unsigned short> offset = Defaults.Offset, bool wrap = Defaults.Wrap, bool padding = Defaults.UsePadding) {return wintp(pos, input, minWidth, color, att, offset, wrap, padding).second;}
            /// @brief Write Integer, Return X-Position - Write an integer to the window - short short pos, short return
            /// @param y y-position (row) of the start of the integer
            /// @param x x-position (col) of the start of the integer
            /// @param input Integer input to be written
            /// @param minWidth Minimum width for the integer to be (leading zeros will be added to meet the width)
            /// @param color Color pair to use when writing each cell
            /// @param att Set of attributes to apply (in any order): bo = Bold, it = Italic, un = Underline, re = Reverse, bl = Blink, di = Dim, in = Invisible, st = Standout, pr = Protected, al = Altset
            /// @param offset Pair consisiting of a y-offset (rows) and an x-offset (cols) that changes the return value
            /// @param wrap Whether to wrap to starting x or not
            /// @param usePadding Whether to disallow writing to cells that are in a padded region or not
            /// @returns A x-position (col) based off of the last character of the integer that was written and the inputted offset
            unsigned short wintx(unsigned short y, unsigned short x, int input, unsigned short minWidth = Defaults.MinWidth, unsigned char color = Defaults.Color, std::string att = Defaults.Attributes, std::pair<unsigned short, unsigned short> offset = Defaults.Offset, bool wrap = Defaults.Wrap, bool padding = Defaults.UsePadding) {return wintx({y, x}, input, minWidth, color, att, offset, wrap, padding);}
            /// @brief Write Integer, Return Nothing - Write an integer to the window - pair pos, no return
            /// @param pos Pair consisting of a y-position (row) and an x-position (col) indicating the start of the integer
            /// @param input Integer input to be written
            /// @param minWidth Minimum width for the integer to be (leading zeros will be added to meet the width)
            /// @param color Color pair to use when writing each cell
            /// @param att Set of attributes to apply (in any order): bo = Bold, it = Italic, un = Underline, re = Reverse, bl = Blink, di = Dim, in = Invisible, st = Standout, pr = Protected, al = Altset
            /// @param wrap Whether to wrap to starting x or not
            /// @param usePadding Whether to disallow writing to cells that are in a padded region or not
            void wint(std::pair<unsigned short, unsigned short> pos, int input, unsigned short minWidth = Defaults.MinWidth, unsigned char color = Defaults.Color, std::string att = Defaults.Attributes, bool wrap = Defaults.Wrap, bool padding = Defaults.UsePadding) {wintp(pos, input, minWidth, color, att, {0, 0}, wrap, padding);}
            /// @brief Write Integer, Return Nothing - Write an integer to the window - short short pos, no return
            /// @param y y-position (row) of the start of the integer
            /// @param x x-position (col) of the start of the integer
            /// @param input Integer input to be written
            /// @param minWidth Minimum width for the integer to be (leading zeros will be added to meet the width)
            /// @param color Color pair to use when writing each cell
            /// @param att Set of attributes to apply (in any order): bo = Bold, it = Italic, un = Underline, re = Reverse, bl = Blink, di = Dim, in = Invisible, st = Standout, pr = Protected, al = Altset
            /// @param wrap Whether to wrap to starting x or not
            /// @param usePadding Whether to disallow writing to cells that are in a padded region or not
            void wint(unsigned short y, unsigned short x, int input, unsigned short minWidth = Defaults.MinWidth, unsigned char color = Defaults.Color, std::string att = Defaults.Attributes, bool wrap = Defaults.Wrap, bool padding = Defaults.UsePadding) {wint({y, x}, input, minWidth, color, att, wrap, padding);}

            /// @brief Write Matrix String, Return Position - Write a string to the window using the 6x6 matrix text - pair pos, pair return
            /// @param pos Pair consisting of a y-position (row) and an x-position (col) indicating the top-left corner of the first character of the string
            /// @param input String input to be written (no unicode, somewhat limited character set)
            /// @param color Color pair to use when writing each cell
            /// @param danglingLetters Whether to allow certain letters to hang beneath the usual line (increases leading by 1)
            /// @param kerning The amount of spaces between each letter
            /// @param leading The amount of spaces between each line
            /// @param offset Pair consisiting of a y-offset (rows) and an x-offset (cols) that changes the return value
            /// @param wrap Whether to wrap to starting x or not
            /// @param usePadding Whether to disallow writing to cells that are in a padded region or not
            /// @returns A pair consisting of a y-position (row) and x-position (col) of the bottom-right corner of the last character of the text that was written
            std::pair<unsigned short, unsigned short> wmstrp(std::pair<unsigned short, unsigned short> pos, std::string input, unsigned char matrix = Defaults.Matrix, unsigned char color = Defaults.Color, bool danglingLetters = Defaults.DanglingLetters, unsigned char kerning = Defaults.Kerning, unsigned char leading = Defaults.Leading, std::pair<unsigned short, unsigned short> offset = Defaults.Offset, bool wrap = Defaults.Wrap, bool padding = Defaults.UsePadding) {
                if (!checkCoord(pos, false)) {return pos;}
                kerning = kerning < 1 ? 1 : kerning;
                leading = leading < 1 ? 1 : leading;
                matrix = matrix < 0 ? 0 : matrix;
                matrix = matrix > Keys.DanglingOffsets.size() - 1 ? Keys.DanglingOffsets.size() - 1 : matrix;

                unsigned short posy = pos.first, posx = 0, xoffset = 0;
                bool dangling = false;
                unsigned char chary = Keys.MatrixText[0].second[matrix].size(), charx = Keys.MatrixText[0].second[matrix][0].length(), danglingOffset = Keys.DanglingOffsets[matrix];

                for (unsigned short i = 0; i < input.length(); i++) {
                    // Automatically wrap the string back to the starting x-pos
                    if (posx + charx > DimX - charx - kerning - (padding ? PadLeft + PadRight : 0)) {
                        posy += chary + leading + (danglingLetters && leading < danglingOffset + 1 ? danglingOffset : 0);
                        xoffset += posx + charx + kerning - (wrap ? pos.second : 0);
                    }

                    posx = pos.second + i * (charx + kerning) - xoffset;

                    // Quit early if the string goes out of bounds
                    if (!checkCoord(posy + chary - 1, posx + charx - 1, false)) {break;}

                    for (unsigned char j = 0; j < Keys.MatrixText.size(); j++) {
                        if (input[i] == Keys.MatrixText[j].first) {
                            if (input[i] == 'g' || input[i] == 'j' || input[i] == 'p' || input[i] == 'q' || input[i] == 'y') {dangling = true;}

                            for (unsigned char k = 0; k < chary; k++) {
                                for (unsigned char l = 0; l < charx; l++) {
                                    wchar(posy + k + (dangling ? danglingOffset : 0), posx + l, Keys.MatrixText[j].second[matrix][k][l], color);
                                }
                            }

                            // Lowercase j needs a bit of extra work since instead of just moving down when dangling it has to get longer
                            if (input[i] == 'j') {
                                for (unsigned char k = 0; k < chary - danglingOffset; k++) {
                                    for (unsigned char l = 0; l < charx; l++) {
                                        wchar(posy + k, posx + l, Keys.MatrixText[j].second[matrix][k][l], color);
                                    }
                                }
                            }

                            dangling = false;
                            break;
                        }
                    }
                }

                return {posy + chary - 1 + offset.first, posx + charx - 1 + offset.second};
            }
            /// @brief Write Matrix String, Return Position - Write a string to the window using the 6x6 matrix text - short short pos, pair return
            /// @param y y-position (row) indicating the top-left corner of the first character of the string
            /// @param x x-position (col) indicating the top-left corner of the first character of the string
            /// @param input String input to be written (no unicode, somewhat limited character set)
            /// @param color Color pair to use when writing each cell
            /// @param danglingLetters Whether to allow certain letters to hang beneath the usual line (increases leading by 1)
            /// @param kerning The amount of spaces between each letter
            /// @param leading The amount of spaces between each line
            /// @param offset Pair consisiting of a y-offset (rows) and an x-offset (cols) that changes the return value
            /// @param wrap Whether to wrap to starting x or not
            /// @param usePadding Whether to disallow writing to cells that are in a padded region or not
            /// @returns A pair consisting of a y-position (row) and x-position (col) of the bottom-right corner of the last character of the text that was written
            std::pair<unsigned short, unsigned short> wmstrp(unsigned short y, unsigned short x, std::string input, unsigned char matrix = Defaults.Matrix, unsigned char color = Defaults.Color, bool danglingLetters = Defaults.DanglingLetters, unsigned char kerning = Defaults.Kerning, unsigned char leading = Defaults.Leading, std::pair<unsigned short, unsigned short> offset = Defaults.Offset, bool wrap = Defaults.Wrap) {return wmstrp({y, x}, input, matrix, color, danglingLetters, kerning, leading, offset, wrap);}
            /// @brief Write Matrix String, Return Y-Position - Write a string to the window using the 6x6 matrix text - pair pos, short return
            /// @param pos Pair consisting of a y-position (row) and an x-position (col) indicating the top-left corner of the first character of the string
            /// @param input String input to be written (no unicode, somewhat limited character set)
            /// @param color Color pair to use when writing each cell
            /// @param danglingLetters Whether to allow certain letters to hang beneath the usual line (increases leading by 1)
            /// @param kerning The amount of spaces between each letter
            /// @param leading The amount of spaces between each line
            /// @param offset Pair consisiting of a y-offset (rows) and an x-offset (cols) that changes the return value
            /// @param wrap Whether to wrap to starting x or not
            /// @param usePadding Whether to disallow writing to cells that are in a padded region or not
            /// @returns The y-position (row) of the bottom-right corner of the last character of the text that was written
            unsigned short wmstry(std::pair<unsigned short, unsigned short> pos, std::string input, unsigned char matrix = Defaults.Matrix, unsigned char color = Defaults.Color, bool danglingLetters = Defaults.DanglingLetters, unsigned char kerning = Defaults.Kerning, unsigned char leading = Defaults.Leading, std::pair<unsigned short, unsigned short> offset = Defaults.Offset, bool wrap = Defaults.Wrap) {return wmstrp(pos, input, matrix, color, danglingLetters, kerning, leading, offset, wrap).first;}
            /// @brief Write Matrix String, Return Y-Position - Write a string to the window using the 6x6 matrix text - short short pos, short return
            /// @param y y-position (row) indicating the top-left corner of the first character of the string
            /// @param x x-position (col) indicating the top-left corner of the first character of the string
            /// @param input String input to be written (no unicode, somewhat limited character set)
            /// @param color Color pair to use when writing each cell
            /// @param danglingLetters Whether to allow certain letters to hang beneath the usual line (increases leading by 1)
            /// @param kerning The amount of spaces between each letter
            /// @param leading The amount of spaces between each line
            /// @param offset Pair consisiting of a y-offset (rows) and an x-offset (cols) that changes the return value
            /// @param wrap Whether to wrap to starting x or not
            /// @param usePadding Whether to disallow writing to cells that are in a padded region or not
            /// @returns The y-position (row) of the bottom-right corner of the last character of the text that was written
            unsigned short wmstry(unsigned short y, unsigned short x, std::string input, unsigned char matrix = Defaults.Matrix, unsigned char color = Defaults.Color, bool danglingLetters = Defaults.DanglingLetters, unsigned char kerning = Defaults.Kerning, unsigned char leading = Defaults.Leading, std::pair<unsigned short, unsigned short> offset = Defaults.Offset, bool wrap = Defaults.Wrap) {return wmstry({y, x}, input, matrix, color, danglingLetters, kerning, leading, offset, wrap);}
            /// @brief Write Matrix String, Return X-Position - Write a string to the window using the 6x6 matrix text - pair pos, short return
            /// @param pos Pair consisting of a y-position (row) and an x-position (col) indicating the top-left corner of the first character of the string
            /// @param input String input to be written (no unicode, somewhat limited character set)
            /// @param color Color pair to use when writing each cell
            /// @param danglingLetters Whether to allow certain letters to hang beneath the usual line (increases leading by 1)
            /// @param kerning The amount of spaces between each letter
            /// @param leading The amount of spaces between each line
            /// @param offset Pair consisiting of a y-offset (rows) and an x-offset (cols) that changes the return value
            /// @param wrap Whether to wrap to starting x or not
            /// @param usePadding Whether to disallow writing to cells that are in a padded region or not
            /// @returns The x-position (col) of the bottom-right corner of the last character of the text that was written
            unsigned short wmstrx(std::pair<unsigned short, unsigned short> pos, std::string input, unsigned char matrix = Defaults.Matrix, unsigned char color = Defaults.Color, bool danglingLetters = Defaults.DanglingLetters, unsigned char kerning = Defaults.Kerning, unsigned char leading = Defaults.Leading, std::pair<unsigned short, unsigned short> offset = Defaults.Offset, bool wrap = Defaults.Wrap) {return wmstrp(pos, input, matrix, color, danglingLetters, kerning, leading, offset, wrap).second;}
            /// @brief Write Matrix String, Return X-Position - Write a string to the window using the 6x6 matrix text - short short pos, short return
            /// @param y y-position (row) indicating the top-left corner of the first character of the string
            /// @param x x-position (col) indicating the top-left corner of the first character of the string
            /// @param input String input to be written (no unicode, somewhat limited character set)
            /// @param color Color pair to use when writing each cell
            /// @param danglingLetters Whether to allow certain letters to hang beneath the usual line (increases leading by 1)
            /// @param kerning The amount of spaces between each letter
            /// @param leading The amount of spaces between each line
            /// @param offset Pair consisiting of a y-offset (rows) and an x-offset (cols) that changes the return value
            /// @param wrap Whether to wrap to starting x or not
            /// @param usePadding Whether to disallow writing to cells that are in a padded region or not
            /// @returns The x-position (col) of the bottom-right corner of the last character of the text that was written
            unsigned short wmstrx(unsigned short y, unsigned short x, std::string input, unsigned char matrix = Defaults.Matrix, unsigned char color = Defaults.Color, bool danglingLetters = Defaults.DanglingLetters, unsigned char kerning = Defaults.Kerning, unsigned char leading = Defaults.Leading, std::pair<unsigned short, unsigned short> offset = Defaults.Offset, bool wrap = Defaults.Wrap) {return wmstry({y, x}, input, matrix, color, danglingLetters, kerning, leading, offset, wrap);}
            /// @brief Write Matrix String, Return Nothing - Write a string to the window using the 6x6 matrix text - pair pos, no return
            /// @param pos Pair consisting of a y-position (row) and an x-position (col) indicating the top-left corner of the first character of the string
            /// @param input String input to be written (no unicode, somewhat limited character set)
            /// @param color Color pair to use when writing each cell
            /// @param danglingLetters Whether to allow certain letters to hang beneath the usual line (increases leading by 1)
            /// @param kerning The amount of spaces between each letter
            /// @param leading The amount of spaces between each line
            /// @param offset Pair consisiting of a y-offset (rows) and an x-offset (cols) that changes the return value
            /// @param wrap Whether to wrap to starting x or not
            /// @param usePadding Whether to disallow writing to cells that are in a padded region or not
            void wmstr(std::pair<unsigned short, unsigned short> pos, std::string input, unsigned char matrix = Defaults.Matrix, unsigned char color = Defaults.Color, bool danglingLetters = Defaults.DanglingLetters, unsigned char kerning = Defaults.Kerning, unsigned char leading = Defaults.Leading, bool wrap = Defaults.Wrap) {wmstrp(pos, input, matrix, color, danglingLetters, kerning, leading, {0, 0}, wrap);}
            /// @brief Write Matrix String, Return Nothing - Write a string to the window using the 6x6 matrix text - short short pos, no return
            /// @param y y-position (row) indicating the top-left corner of the first character of the string
            /// @param x x-position (col) indicating the top-left corner of hte first character of the string
            /// @param input String input to be written (no unicode, somewhat limited character set)
            /// @param color Color pair to use when writing each cell
            /// @param danglingLetters Whether to allow certain letters to hang beneath the usual line (increases leading by 1)
            /// @param kerning The amount of spaces between each letter
            /// @param leading The amount of spaces between each line
            /// @param offset Pair consisiting of a y-offset (rows) and an x-offset (cols) that changes the return value
            /// @param wrap Whether to wrap to starting x or not
            /// @param usePadding Whether to disallow writing to cells that are in a padded region or not
            void wmstr(unsigned short y, unsigned short x, std::string input, unsigned char matrix = Defaults.Matrix, unsigned char color = Defaults.Color, bool danglingLetters = Defaults.DanglingLetters, unsigned char kerning = Defaults.Kerning, unsigned char leading = Defaults.Leading, bool wrap = Defaults.Wrap) {wmstr({y, x}, input, matrix, color, danglingLetters, kerning, leading, wrap);}

            //
            // GET USER INPUT
            //

            /// @brief Get Character - Get a single character input from the user (acts as the ncursespp version of wgetch())
            /// @param enableKeypad Allow the use of arrow keys and such
            /// @param pause Pause the program until an input is read
            /// @param autoRender Automatically render the window when the function is called
            /// @returns An integer that can be used in the same way that it is in the ncurses getch() function
            int gchar(bool pause = true, bool enableKeypad = true, bool autoRender = true) {
                if (autoRender) {rinst();}

                if (enableKeypad) {keypad(Win, true);}
                else {keypad(Win, false);}
                if (pause) {nodelay(Win, false);}
                else {nodelay(Win, true);}

                return wgetch(Win);
            }

            /// @brief Get String - Get a string input from the user (acts as the ncursespp version of mvwgetnstr())
            /// @param y y-position (row) of the start of the input location
            /// @param x x-position (col) of the start of the input location
            /// @param maxChars Maximum amount of characters the string can contain
            /// @param echoColor Color pair to use when echoing the user's input
            /// @param echoAtt Set of attributes to apply to the echoed string
            /// @param autoWrite Whether to automatically write the inputted string to the window
            /// @param showStr Whether to echo the user's input or not
            /// @param showCursor Whether to show the cursor when taking an input or not
            /// @param enableKeypad Whether to allow the use of the keypad when taking the input
            /// @returns A wide string (not unicode friendly :/)
            std::wstring gstr(unsigned short y, unsigned short x, int maxChars = 255, unsigned char echoColor = Defaults.Color, std::string echoAtt = Defaults.Attributes, bool autoWrite = true, bool showStr = true, bool showCursor = true, bool enableKeypad = true) {
                wint_t input[maxChars];
                std::wstring output;

                if (enableKeypad) {keypad(Win, true);}
                if (showStr) {echo();}
                if (showCursor) {curs_set(1);}

                wattron(Win, COLOR_PAIR(echoColor));
                toggleAttributes(echoAtt);

                mvwgetn_wstr(Win, y, x, input, maxChars);

                wattroff(Win, COLOR_PAIR(echoColor));
                toggleAttributes(echoAtt);

                for (int i = 0; i < maxChars; i++) {
                    if (input[i] == '\0') {break;}
                    output += input[i];
                }

                if (autoWrite) {wstr(y, x, output, echoColor, echoAtt);}

                curs_set(0);
                noecho();
                keypad(Win, false);

                return output;
            }

            //
            // RENDERING THE WINDOW
            //

            /// @brief Render Instantly - Render the window instantly
            void rinst() {
                for (unsigned short i = 0; i < DimY; i++) {
                    for (unsigned short j = 0; j < DimX; j++) {
                        write(i, j);
                    }
                }
                wrefresh(Win);
            }
            /// @brief Render Instantly - Render a portion of the window instantly (mainly for performance on weaker computers)
            /// @param y y-position (row) to start rendering the window at
            /// @param x x-position (col) to start rendering the window at
            /// @param dimy Amount of rows to render
            /// @param dimx Amount of columns to render
            /// @param usePadding Whether to offset coordinates and disallow the editing of certain cells or not
            void rinst(unsigned short y, unsigned short x, unsigned short dimy, unsigned short dimx, bool usePadding = Defaults.UsePadding) {
                if (!checkCoord(y, x, usePadding) || !checkCoord(y + dimy - 1, x + dimx - 1, usePadding)) {return;}
                
                for (unsigned short i = 0; i < dimy; i++) {
                    for (unsigned short j = 0; j < dimx; j++) {
                        write(y + i, x + j);
                    }
                }
                wrefresh(Win);
            }

            /// @brief Render by Line - Render the window line-by-line (and char-by-char if indicated)
            /// @param dir Direction to start from: 0 = Top-Left, 1 = Bottom-Left, 2 = Left-Top, 3 = Right-Top
            /// @param full Whether to render the full line at once or not
            /// @param rev Whether to reverse the char-by-char direction of rendering (has no effect if full is true)
            /// @param millis Milliseconds to wait in between each line/character rendering
            void rline(unsigned char dir = DIR_UP, bool full = true, bool rev = false, bool ignorePadding = Defaults.UsePadding, unsigned long millis = 20) {
                if (dir < 0 || dir > 3) {return rinst();}

                unsigned short l1 = dir < 2 ? (DimY - (ignorePadding ? (rev ? PadUp : PadDown) : 0)) : (DimX - (ignorePadding ? (rev ? PadLeft : PadRight) : 0));
                unsigned short l2 = dir < 2 ? (DimX - (ignorePadding ? (rev ? PadLeft : PadRight) : 0)) : (DimY - (ignorePadding ? (rev ? PadUp : PadDown) : 0));
                unsigned short o1 = dir < 2 ? (rev ? (ignorePadding ? PadDown : 0) : (ignorePadding ? PadUp : 0)) : (rev ? (ignorePadding ? PadRight : 0) : (ignorePadding ? PadLeft : 0));
                unsigned short o2 = dir < 2 ? (rev ? (ignorePadding ? PadRight : 0) : (ignorePadding ? PadLeft : 0)) : (rev ? (ignorePadding ? PadDown : 0) : (ignorePadding ? PadUp : 0));
                
                for (unsigned short i = o1; i < l1; i++) {
                    for (unsigned short j = o2; j < l2; j++) {
                        switch (dir) {
                            case 0:
                                write(i, rev ? DimX - 1 - j : j);
                                break;
                            case 1:
                                write(DimY - 1 - i, rev ? DimX - 1 - j : j);
                                break;
                            case 2:
                                write(rev ? DimY - 1 - j : j, i);
                                break;
                            case 3:
                                write(rev ? DimY - 1 - j : j, DimX - 1 -i);
                                break;
                        }

                        if (!full) {
                            wrefresh(Win);
                            if (wait(millis)) {return rinst();}
                        }
                    }
                    
                    wrefresh(Win);
                    if (wait(millis)) {return rinst();}
                }
            }
            /// @brief Render by Line (from Top) - Render the window line-by-line (and char-by-char if indicated) from the top of the window
            /// @param full Whether to render the full line at once or not
            /// @param rev Whether to reverse the char-by-char direction of rendering (has no effect if full is true)
            /// @param millis Milliseconds to wait in between each line/character rendering
            void rlinetop(bool full = true, bool rev = false, bool padding = Defaults.UsePadding, unsigned long millis = 20) {rline(0, full, rev, padding, millis);}
            /// @brief Render by Line (from Bottom) - Render the window line-by-line (and char-by-char if indicated) from the bottom of the window
            /// @param full Whether to render the full line at once or not
            /// @param rev Whether to reverse the char-by-char direction of rendering (has no effect if full is true)
            /// @param millis Milliseconds to wait in between each line/character rendering
            void rlinebot(bool full = true, bool rev = false, bool padding = Defaults.UsePadding, unsigned long millis = 20) {rline(1, full, rev, padding, millis);}
            /// @brief Render by Line (from Left) - Render the window line-by-line (and char-by-char if indicated) from the left of the window
            /// @param full Whether to render the full line at once or not
            /// @param rev Whether to reverse the char-by-char direction of rendering (has no effect if full is true)
            /// @param millis Milliseconds to wait in between each line/character rendering
            void rlineleft(bool full = true, bool rev = false, bool padding = Defaults.UsePadding, unsigned long millis = 20) {rline(2, full, rev, padding, millis);}
            /// @brief Render by Line (from Right) - Render the window line-by-line (and char-by-char if indicated) from the right of the window
            /// @param full Whether to render the full line at once or not
            /// @param rev Whether to reverse the char-by-char direction of rendering (has no effect if full is true)
            /// @param millis Milliseconds to wait in between each line/character rendering
            void rlineright(bool full = true, bool rev = false, bool padding = Defaults.UsePadding, unsigned long millis = 20) {rline(3, full, rev, padding, millis);}

            /// @brief Render in a Radial Motion - Render the window by "sweeping" a line in a radial motion and rendering the window as the line passes over each cell
            /// @param divisions Amount of divisions to split the window into and render individually
            /// @param angle The starting angle of the main "sweeper" (expressed in degrees, 0 = midpoint of the right edge of the window)
            /// @param ccw Whether to "sweep" in a counter-clockwise direction or not
            /// @param millis Milliseconds to wait in between each line/character rendering
            /// @param resolution Number that basically dictates how many steps are required to render each division (smaller number = slower, more accurate "sweeps")
            void rrad(unsigned char divisions = 2, double angle = 90, bool ccw = true, unsigned long millis = 5, double resolution = 0.005) {
                // Fix the angle input and then convert it to radians
                angle = angle < 0 && !ccw ? angle + 90 : angle;
                while (angle < 0) {angle += 360;}
                angle = angle * (M_PI / 180);

                double divAngle = 2 * M_PI / divisions;
                double slope, cangle;
                unsigned char resMult;

                for (double i = angle; i < divAngle + angle; i += resolution) {
                    // Each division will need to draw its own line and stuff
                    for (unsigned char j = 0; j < divisions; j++) {
                        // Find the angle that will be used to find the current slope
                        cangle = i + divAngle * j;
                        while (cangle >= M_PI * 2) {cangle -= M_PI * 2;}

                        slope = (ccw ? -1 : 1) * ((DimX / 2) * sin(cangle)) / (DimX * cos(cangle));
                        resMult = (cangle < M_PI / 2 + resolution * 16 && cangle > M_PI / 2 - resolution * 16) || (cangle < 3 * M_PI / 2 + resolution * 16 && cangle > 3 * M_PI / 2 - resolution * 16) ? 10 : 25;

                        // Pretty much just draw a line from the middle of the window to the edge of the window along the calculated slope
                        for (double k = 0; k <= DimX / 2; k+= resolution * resMult) {
                            if (cangle < M_PI / 2 || cangle > 3 * M_PI / 2) {
                                write(DimY / 2 + slope * k, DimX / 2 + k * 2);
                                write(DimY / 2 + slope * k - 1, DimX / 2 + k * 2);
                                write(DimY / 2 + slope * k + 1, DimX / 2 + k * 2);
                                write(DimY / 2 + slope * k, DimX / 2 + k * 2 - 1);
                                write(DimY / 2 + slope * k, DimX / 2 + k * 2 + 2);
                            } else {
                                write(DimY / 2 - slope * k, DimX / 2 - k * 2);
                                write(DimY / 2 - slope * k - 1, DimX / 2 - k * 2);
                                write(DimY / 2 - slope * k + 1, DimX / 2 - k * 2);
                                write(DimY / 2 - slope * k, DimX / 2 - k * 2 - 1);
                                write(DimY / 2 - slope * k, DimX / 2 - k * 2 + 2);
                            }
                        }
                    }
                
                    wrefresh(Win);
                    if (wait(millis)) {return rinst();}
                }

                rinst();
            }

            //
            // LINE DRAWING
            //

            /// @brief Draw Line - Draw a line with unicode box drawing characters
            /// @param y Starting y-position (row) for the line
            /// @param x Starting x-position (col) for the line
            /// @param length Length of the line
            /// @param vertical Whether the line is vertical or not
            /// @param rev Whether the line is reversed or not (normal is up to down or left to right and reverse switches the direction)
            /// @param style A pair consisting of a number corresponding to the main style (0 = Light Hard, 1 = Light Soft, 2 = Heavy Both, 3 = Heavy Horizontal, 4 = Heavy Vertical, 5 = Doubled Both, 6 = Doubled Horizontal, 7 = Doubled Vertical) and a number corresponding to the dash style (0 = None, 1 = Double, 2 = Triple, 3 = Quardruple, 4 = Single 1, 5 = Single 2)
            /// @param color Color pair to use when drawing the line
            /// @param mergeable Whether the line will be able to be merged with future lines
            /// @param canMerge Whether the line will be able to merge with pre-existing box-drawing characters
            void dline(unsigned short y, unsigned short x, unsigned short length, bool vertical = false, bool rev = false, std::pair<unsigned char, unsigned char> style = Defaults.Style, unsigned char color = Defaults.Color, bool padding = Defaults.UsePadding, bool mergeable = Defaults.Mergeable, bool canMerge = Defaults.CanMerge) {
                if (!checkCoord(y, x, false) || length <= 0 || (y < 0 && rev) || (y >= DimY && !rev) || (x < 0 && rev) || (x >= DimX && !rev)) {return;}
                style.first = (style.first < 0 || style.first >= 8) ? LIGHT_HARD : style.first;
                style.second = (style.second < 0 || style.second >= 6) ? DASHED_NONE : style.second;

                unsigned short posy = y, posx = x;
                std::vector<unsigned char> dirs;
                wchar_t piece;
                unsigned char tdir1 = vertical ? (DIR_DOWN - (rev ? 1 : 0)) : (DIR_RIGHT - (rev ? 1 : 0));
                unsigned char tdir2 = tdir1 + (rev ? 1 : -1);

                /// @brief Used to fabricate a direction (connection style) based on the inputted style
                unsigned char cstyle;
                switch (style.first) {
                    case LIGHT_HARD:
                    case LIGHT_SOFT:
                        cstyle = STYLE_LIGHT;
                        break;
                    case HEAVY_BOTH:
                        cstyle = STYLE_HEAVY;
                        break;
                    case HEAVY_HORIZONTAL:
                        cstyle = vertical ? STYLE_LIGHT : STYLE_HEAVY;
                        break;
                    case HEAVY_VERTICAL:
                        cstyle = vertical ? STYLE_HEAVY : STYLE_LIGHT;
                        break;
                    case DOUBLED_BOTH:
                        cstyle = STYLE_DOUBLED;
                        break;
                    case DOUBLED_HORIZTONAL:
                        cstyle = vertical ? STYLE_LIGHT : STYLE_DOUBLED;
                        break;
                    case DOUBLED_VERTICAL:
                        cstyle = vertical ? STYLE_DOUBLED : STYLE_LIGHT;
                        break;
                }

                for (unsigned short i = 0; i < length; i++) {
                    // Quit early if the line starts to venture out of the window
                    if (!checkCoord(posy, posx, false)) {break;}

                    // Get the directions (connection styles) of the surrounding characters (for connecting to existing lines)
                    // This step is skipped if the line isn't suppsoed to merge with others
                    if (canMerge) {dirs = {getConnectStyle(posy - 1, posx, DIR_DOWN, padding), getConnectStyle(posy + 1, posx, DIR_UP, padding), getConnectStyle(posy, posx - 1, DIR_RIGHT, padding), getConnectStyle(posy, posx + 1, DIR_LEFT, padding)};}
                    else {dirs = {STYLE_NONE, STYLE_NONE, STYLE_NONE, STYLE_NONE};}

                    // With the exception of the end pieces, the line has to have directions fabricated in order to actually exist
                    if (i == 0) {dirs[tdir2] = (dirs[tdir2] == STYLE_NONE && dirs[vertical ? DIR_LEFT : DIR_UP] == STYLE_NONE && dirs[vertical ? DIR_RIGHT : DIR_DOWN] == STYLE_NONE) ? cstyle : dirs[tdir2];}
                    else {dirs[tdir2] = cstyle;}
                    if (i == length - 1) {dirs[tdir1] = (dirs[tdir1] == STYLE_NONE && dirs[vertical ? DIR_LEFT : DIR_UP] == STYLE_NONE && dirs[vertical ? DIR_RIGHT : DIR_DOWN] == STYLE_NONE) ? cstyle : dirs[tdir1];}
                    else {dirs[tdir1] = cstyle;}

                    piece = getPiece(dirs, style);

                    // In the scenario that a doubled line intersects with a heavy one, the directions need to be altered because there are no
                    // unicode characters that properly mesh the two styles; the current style is kept because it is the most recent one being used
                    if (piece == L' ') {
                        // If the current style is heavy, then all doubled directions are replaced with heavy directions
                        // The opposite is done if the current style is doubled
                        for (unsigned char i = 0; i < 4; i++) {
                            if (style.first == HEAVY_BOTH || style.first == HEAVY_HORIZONTAL || style.first == HEAVY_VERTICAL) {dirs[i] = (dirs[i] == STYLE_DOUBLED) ? STYLE_HEAVY : dirs[i];}
                            else {dirs[i] = (dirs[i] == STYLE_HEAVY) ? STYLE_DOUBLED : dirs[i];}
                        }

                        piece = getPiece(dirs, style);
                    }

                    wchar(posy, posx, piece, color, "", padding);
                    Grid[posy + (padding ? PadUp : 0)][posx + (padding ? PadLeft : 0)].CanMerge = mergeable;

                    posy = y + (rev ? -1 : 1) * (i + 1) * (vertical ? 1 : 0);
                    posx = x + (rev ? -1 : 1) * (i + 1) * (vertical ? 0 : 1);
                }
            }
            /// @brief Draw Vertical Line - Draw a vertical line (automatically merges with pre-existing lines)
            /// @param y Starting y-position (row) for the line
            /// @param x Starting x-position (col) for the line
            /// @param length Length of the line
            /// @param rev Whether the line is reversed or not (normal is up to down or left to right and reverse switches the direction)
            /// @param style A pair consisting of a number corresponding to the main style (0 = Light Hard, 1 = Light Soft, 2 = Heavy Both, 3 = Heavy Horizontal, 4 = Heavy Vertical, 5 = Doubled Both, 6 = Doubled Horizontal, 7 = Doubled Vertical) and a number corresponding to the dash style (0 = None, 1 = Double, 2 = Triple, 3 = Quardruple, 4 = Single 1, 5 = Single 2)
            /// @param color Color pair to use when drawing the line
            /// @param mergeable Whether the line will be able to be merged with future lines
            /// @param canMerge Whether the line will be able to merge with pre-existing box-drawing characters
            void dvline(unsigned short y, unsigned short x, unsigned short length, bool rev = false, std::pair<unsigned char, unsigned char> style = Defaults.Style, unsigned char color = Defaults.Color, bool padding = Defaults.UsePadding, bool mergeable = Defaults.Mergeable, bool canMerge = Defaults.CanMerge) {dline(y, x, length, true, rev, style, color, padding, mergeable, canMerge);}
            /// @brief Draw Horizontal Line - Draw a horizontal line (automatically merges with pre-existing lines)
            /// @param y Starting y-position (row) for the line
            /// @param x Starting x-position (col) for the line
            /// @param length Length of the line
            /// @param rev Whether the line is reversed or not (normal is up to down or left to right and reverse switches the direction)
            /// @param style A pair consisting of a number corresponding to the main style (0 = Light Hard, 1 = Light Soft, 2 = Heavy Both, 3 = Heavy Horizontal, 4 = Heavy Vertical, 5 = Doubled Both, 6 = Doubled Horizontal, 7 = Doubled Vertical) and a number corresponding to the dash style (0 = None, 1 = Double, 2 = Triple, 3 = Quardruple, 4 = Single 1, 5 = Single 2)
            /// @param color Color pair to use when drawing the line
            /// @param mergeable Whether the line will be able to be merged with future lines
            /// @param canMerge Whether the line will be able to merge with pre-existing box-drawing characters
            void dhline(unsigned short y, unsigned short x, unsigned short length, bool rev = false, std::pair<unsigned char, unsigned char> style = Defaults.Style, unsigned char color = Defaults.Color, bool padding = Defaults.UsePadding, bool mergeable = Defaults.Mergeable, bool canMerge = Defaults.CanMerge) {dline(y, x, length, false, rev, style, color, padding, mergeable, canMerge);}

            /// @brief Draw Box - Draw a box - top left
            /// @param y y-position (row) of the top-left corner of the box
            /// @param x x-position (col) of the top-left corner of the box
            /// @param dimy Height (rows) of the box
            /// @param dimx Length (cols) of the box
            /// @param style A pair consisting of a number corresponding to the main style (0 = Light Hard, 1 = Light Soft, 2 = Heavy Both, 3 = Heavy Horizontal, 4 = Heavy Vertical, 5 = Doubled Both, 6 = Doubled Horizontal, 7 = Doubled Vertical) and a number corresponding to the dash style (0 = None, 1 = Double, 2 = Triple, 3 = Quardruple, 4 = Single 1, 5 = Single 2)
            /// @param color Color pair to use when drawing the line
            /// @param mergeable Whether the line will be able to be merged with future lines
            /// @param canMerge Whether the line will be able to merge with pre-existing box-drawing characters
            void dbox(unsigned short y, unsigned short x, unsigned short dimy, unsigned short dimx, std::pair<unsigned char, unsigned char> style = Defaults.Style, unsigned char color = Defaults.Color, bool padding = Defaults.UsePadding, bool mergeable = Defaults.Mergeable, bool canMerge = Defaults.CanMerge) {                
                if (!checkCoord(y, x, false) || !checkCoord(y + dimy - 1, x + dimx - 1, false) || dimy < 1 || dimx < 1) {return;}

                // Special cases for when the dimensions of the box are 1
                if (dimy == 1 && dimx == 1) {return wchar(y, x, L'□', color, "", padding);}
                else if (dimy == 1) {return dhline(y, x, dimx, false, style, color, padding, mergeable, canMerge);}
                else if (dimx == 1) {return dvline(y, x, dimy, false, style, color, padding, mergeable, canMerge);}

                dvline(y, x, dimy, false, style, color, padding, mergeable, canMerge);
                dvline(y, x + dimx - 1, dimy, false, style, color, padding, mergeable, canMerge);
                dhline(y, x, dimx, false, style, color, padding, mergeable, canMerge);
                dhline(y + dimy - 1, x, dimx, false, style, color, padding, mergeable, canMerge);
            }
            /// @brief Draw Box - Draw a box - border
            /// @param style A pair consisting of a number corresponding to the main style (0 = Light Hard, 1 = Light Soft, 2 = Heavy Both, 3 = Heavy Horizontal, 4 = Heavy Vertical, 5 = Doubled Both, 6 = Doubled Horizontal, 7 = Doubled Vertical) and a number corresponding to the dash style (0 = None, 1 = Double, 2 = Triple, 3 = Quardruple, 4 = Single 1, 5 = Single 2)
            /// @param color Color pair to use when drawing the line
            /// @param mergeable Whether the line will be able to be merged with future lines
            /// @param canMerge Whether the line will be able to merge with pre-existing box-drawing characters
            void dbox(std::pair<unsigned char, unsigned char> style = Defaults.Style, unsigned char color = Defaults.Color, bool padding = Defaults.UsePadding, bool mergeable = Defaults.Mergeable, bool canMerge = Defaults.CanMerge) {dbox(0, 0, DimY, DimX, style, color, padding, mergeable, canMerge);}
            /// @brief Draw Centered Box - Draw a box - center
            /// @param y y-position (row) of the center of the box
            /// @param x x-position (col) of the center of the box
            /// @param dimy Height (rows) of the box
            /// @param dimx Length (cols) of the box
            /// @param style A pair consisting of a number corresponding to the main style (0 = Light Hard, 1 = Light Soft, 2 = Heavy Both, 3 = Heavy Horizontal, 4 = Heavy Vertical, 5 = Doubled Both, 6 = Doubled Horizontal, 7 = Doubled Vertical) and a number corresponding to the dash style (0 = None, 1 = Double, 2 = Triple, 3 = Quardruple, 4 = Single 1, 5 = Single 2)
            /// @param color Color pair to use when drawing the line
            /// @param mergeable Whether the line will be able to be merged with future lines
            /// @param canMerge Whether the line will be able to merge with pre-existing box-drawing characters
            void dcbox(unsigned short y, unsigned short x, unsigned short dimy, unsigned short dimx, std::pair<unsigned char, unsigned char> style = Defaults.Style, unsigned char color = Defaults.Color, bool padding = Defaults.UsePadding, bool mergeable = Defaults.Mergeable, bool canMerge = Defaults.CanMerge) {dbox(y - dimy / 2, x - dimx / 2, dimy, dimx, style, color, padding, mergeable, canMerge);}

            /// @brief Draw Grid - Draw a grid - top left
            /// @param y y-position of the top-left corner of the grid
            /// @param x x-position of the top-left corner of the grid
            /// @param rows Amount of rows in the grid
            /// @param cols Amount of columns in the grid
            /// @param celly Amount of characters between the horizontal lines that make the cell borders
            /// @param cellx Amount of characters between the vertical lines that make the cell borders
            /// @param style A pair consisting of a number corresponding to the main style (0 = Light Hard, 1 = Light Soft, 2 = Heavy Both, 3 = Heavy Horizontal, 4 = Heavy Vertical, 5 = Doubled Both, 6 = Doubled Horizontal, 7 = Doubled Vertical) and a number corresponding to the dash style (0 = None, 1 = Double, 2 = Triple, 3 = Quardruple, 4 = Single 1, 5 = Single 2)
            /// @param color Color pair to use when drawing the line
            /// @param mergeable Whether the line will be able to be merged with future lines
            /// @param canMerge Whether the line will be able to merge with pre-existing box-drawing characters
            void dgrid(unsigned short y, unsigned short x, unsigned short rows, unsigned short cols, unsigned short celly, unsigned short cellx, std::pair<unsigned char, unsigned char> style = Defaults.Style, unsigned char color = Defaults.Color, bool padding = Defaults.UsePadding, bool mergeable = Defaults.Mergeable, bool canMerge = true) {
                std::pair<unsigned short, unsigned short> dims = gGridDims(rows, cols, celly, cellx);

                if (!checkCoord(y, x, false) || !checkCoord(dims.first, dims.second, false) || rows < 1 || cols < 1 || celly < 0 || cellx < 0) {return;}

                // Vertical lines
                for (unsigned short i = 0; i <= cols; i++) {
                    dvline(y, x + cellx * i + i, dims.first, false, style, color, padding, mergeable, canMerge);
                }
                // Horizontal lines
                for (unsigned short i = 0; i <= rows; i++) {
                    dhline(y + celly * i + i, x, dims.second, false, style, color, padding, mergeable, canMerge);
                }
            }
            /// @brief Draw Grid - Draw a square grid - top left
            /// @param y y-position of the top-left corner of the grid
            /// @param x x-position of the top-left corner of the grid
            /// @param cells Amount of cells in the x and y direction
            /// @param cellDim Amount of characters between the horizontal lines that make the cell borders (this number is doubled for the x-dim of each cell)
            /// @param style A pair consisting of a number corresponding to the main style (0 = Light Hard, 1 = Light Soft, 2 = Heavy Both, 3 = Heavy Horizontal, 4 = Heavy Vertical, 5 = Doubled Both, 6 = Doubled Horizontal, 7 = Doubled Vertical) and a number corresponding to the dash style (0 = None, 1 = Double, 2 = Triple, 3 = Quardruple, 4 = Single 1, 5 = Single 2)
            /// @param color Color pair to use when drawing the line
            /// @param mergeable Whether the line will be able to be merged with future lines
            /// @param canMerge Whether the line will be able to merge with pre-existing box-drawing characters
            void dgrid(unsigned short y, unsigned short x, unsigned short cells, unsigned short cellDim, std::pair<unsigned char, unsigned char> style = Defaults.Style, unsigned char color = Defaults.Color, bool padding = Defaults.UsePadding, bool mergeable = Defaults.Mergeable, bool canMerge = Defaults.CanMerge) {dgrid(y, x, cells, cells, cellDim, cellDim * 2, style, color, padding, mergeable, canMerge);}
            /// @brief Draw Centered Grid - Draw a grid - center
            /// @param y y-position of the center corner of the grid
            /// @param x x-position of the center corner of the grid
            /// @param rows Amount of rows in the grid
            /// @param cols Amount of columns in the grid
            /// @param celly Amount of characters between the horizontal lines that make the cell borders
            /// @param cellx Amount of characters between the vertical lines that make the cell borders
            /// @param style A pair consisting of a number corresponding to the main style (0 = Light Hard, 1 = Light Soft, 2 = Heavy Both, 3 = Heavy Horizontal, 4 = Heavy Vertical, 5 = Doubled Both, 6 = Doubled Horizontal, 7 = Doubled Vertical) and a number corresponding to the dash style (0 = None, 1 = Double, 2 = Triple, 3 = Quardruple, 4 = Single 1, 5 = Single 2)
            /// @param color Color pair to use when drawing the line
            /// @param mergeable Whether the line will be able to be merged with future lines
            /// @param canMerge Whether the line will be able to merge with pre-existing box-drawing characters
            void dcgrid(unsigned short y, unsigned short x, unsigned short rows, unsigned short cols, unsigned short celly, unsigned short cellx, std::pair<unsigned char, unsigned char> style = Defaults.Style, unsigned char color = Defaults.Color, bool padding = Defaults.UsePadding, bool mergeable = Defaults.Mergeable, bool canMerge = Defaults.CanMerge) {dgrid(y - gGridDims(rows, cols, celly, cellx).first, x - gGridDims(rows, cols, celly, cellx).second, rows, cols, celly, cellx, style, color, padding, mergeable, canMerge);}
            /// @brief Draw Centered Grid - Draw a square grid - center
            /// @param y y-position of the center of the grid
            /// @param x x-position of the center of the grid
            /// @param cells Amount of cells in the x and y direction
            /// @param cellDim Amount of characters between the horizontal lines that make the cell borders (this number is doubled for the x-dim of each cell)
            /// @param style A pair consisting of a number corresponding to the main style (0 = Light Hard, 1 = Light Soft, 2 = Heavy Both, 3 = Heavy Horizontal, 4 = Heavy Vertical, 5 = Doubled Both, 6 = Doubled Horizontal, 7 = Doubled Vertical) and a number corresponding to the dash style (0 = None, 1 = Double, 2 = Triple, 3 = Quardruple, 4 = Single 1, 5 = Single 2)
            /// @param color Color pair to use when drawing the line
            /// @param mergeable Whether the line will be able to be merged with future lines
            /// @param canMerge Whether the line will be able to merge with pre-existing box-drawing characters
            void dcgrid(unsigned short y, unsigned short x, unsigned short cells, unsigned short cellDim, std::pair<unsigned char, unsigned char> style = Defaults.Style, unsigned char color = Defaults.Color, bool padding = Defaults.UsePadding, bool mergeable = Defaults.Mergeable, bool canMerge = Defaults.CanMerge) {dcgrid(y, x, cells, cells, cellDim, cellDim * 2, style, color, padding, mergeable, canMerge);}

            /// @brief Get Grid Dimensions - Figure out the x and y dimensions of a grid with the given amount of rows/cols and cell sizes
            /// @param rows Amount of rows in the grid
            /// @param cols Amount of columns in the grid
            /// @param celly Amount of characters between the horizontal lines that make the cell borders
            /// @param cellx Amount of characters between the vertical lines that make the cell borders
            /// @returns A pair consisting of a y-dimension (rows) and x-dimension (cols)
            std::pair<unsigned short, unsigned short> gGridDims(unsigned short rows, unsigned short cols, unsigned short celly, unsigned short cellx) {return {rows * celly + (rows + 1), cols * cellx + (cols + 1)};}
    } mwin;

    class Button {
        private:
            /// @brief Y-dimension (rows) of the button
            unsigned short DimY;
            /// @brief X-dimension (cols) of the button
            unsigned short DimX;
            /// @brief Y-position (row) of the top-left corner of the button
            unsigned short PosY;
            /// @brief X-position (col) of the top-left corner of the button
            unsigned short PosX;

            bool Clickable = true;

            std::vector<char> Targets;

            /// @brief Check Coordinate - Check to see if the mouse position is within the button's bounds
            /// @returns True if the mouse is in bounds, false if it is out of bounds
            bool checkCoord() {return !(Mouse.gy() < PosY || Mouse.gy() >= PosY + DimY || Mouse.gx() < PosX || Mouse.gx() >= PosX + DimX);}

        public:
            Button(unsigned short y = 0, unsigned short x = 0, unsigned short dimy = LINES, unsigned short dimx = COLS, std::vector<char> targets = {}) {
                PosY = (y < 0 || y >= LINES) ? 0 : y;
                PosX = (x < 0 || x >= COLS) ? 0 : x;
                DimY = (dimy < 1 || LINES - dimy - y < 0) ? LINES - y : dimy;
                DimX = (dimx < 1 || COLS - dimx - x < 0) ? COLS - x : dimx;

                tset(targets);
            }
            Button (Window &win, std::vector<char> targets = {}) {
                PosY = win.gposy();
                PosX = win.gposx();
                DimY = win.gdimy();
                DimX = win.gdimx();

                tset(targets);
            }

            void move(unsigned short y, unsigned short x) {
                PosY = (y < 0) ? 0 : y;
                PosY = (y + DimY > LINES) ? LINES - DimY : y;
                PosX = (x < 0) ? 0 : x;
                PosX = (x + DimX > COLS) ? COLS - DimX : x;
            }

            void resize(unsigned short dimy, unsigned short dimx) {
                DimY = (dimy < 1) ? 1 : dimy;
                DimY = (dimy + PosY > LINES) ? LINES - PosY : dimy;
                DimX = (dimx < 1) ? 1 : dimx;
                DimX = (dimx + PosX > COLS) ? COLS - PosX : dimx;
            }

            /// @brief Get Y-Dimension - Get the y-dimension (rows) of a button
            /// @returns The y-dimension (rows) of a button
            const unsigned short gdimy() {return DimY;}
            /// @brief Get X-Dimension - Get the x-dimension (cols) of a button
            /// @returns The x-dimension (cols) of a button
            const unsigned short gdimx() {return DimX;}
            /// @brief Get Y-Position - Get the y-position (row) of the top-left corner of a button
            /// @returns The y-dimension (row) of the top-left corner of a button
            const unsigned short gposy() {return PosY;}
            /// @brief Get X-Position - Get the x-position (col) of the top-left corner of a button
            /// @returns The x-dimension (col) of the top-left corner of a button
            const unsigned short gposx() {return PosX;}
    
            void uclick(bool canClick) {canClick = canClick;}

            /// @brief Check Click - Check to see if the button was pressed or not
            /// @returns The detected mouse event if the button was pressed or M_UNKNOWN (-1) if the button wasn't
            const char cclick() {
                if (!Clickable || !checkCoord() || Mouse.ginput() == M_UNKNOWN) {return M_UNKNOWN;}

                for (unsigned char i = 0; i < Targets.size(); i++) {
                    if (Mouse.ginput() == Targets[i]) {return Targets[i];}
                }

                return M_UNKNOWN;
            }

            /// @brief Target Add - Add a mouse button for the button to look for
            /// @param target A ncursespp mouse button (such as M1_CLICK)
            /// @returns True if the new target is added to the list of targets, false if not
            bool tadd(char target) {
                if (target < 0 || target >= 25) {return false;}

                for (unsigned char i = 0; i < Targets.size(); i++) {
                    if (target == Targets[i]) {return true;}
                }

                Targets.emplace_back(target);
                return true;
            }

            bool tremove(char target) {
                if (target < 0 || target >= 25) {return false;}

                for (unsigned char i = 0; i < Targets.size(); i++) {
                    if (target == Targets[i]) {
                        Targets.erase(Targets.begin() + i);
                        return true;
                    }
                }

                return false;
            }

            bool tset(std::vector<char> targets) {
                std::vector<char> buff;

                for (unsigned char i = 0; i < targets.size(); i++) {
                    if (targets[i] < 0 || targets[i] >= 25) {return false;}
                    buff.emplace_back(targets[i]);
                }

                Targets = buff;
                return true;
            }
    };

    /// @brief Initialize - Acts as the ncursespp version of initscr() with a few other initializations; end() must be called at the end of a program
    /// @param useMouse Whether to set up the ability to take mouse inputs or not
    /// @returns True or false for whether the user's terminal will have full color support
    bool init(bool useMouse = false) {
        setlocale(LC_ALL, "");
        initscr();
        noecho();
        curs_set(0);
        raw();
        srand(time(NULL));

        if (!can_change_color() || !has_colors()) {
            mvwprintw(stdscr, LINES / 2 - 1, COLS / 2 - 50, "Warning: Your terminal doesn't support full color capabilites...");
            mvwprintw(stdscr, LINES / 2, COLS / 2 - 50, "ncursespp requires a terminal that allows for the changing of colors, so go find one that will work");
            mvwprintw(stdscr, LINES / 2 + 1, COLS / 2 - 50, "Press any key to continue...");
            wgetch(stdscr);
            wclear(stdscr);

            return false;
        }

        start_color();

        init_color(NPP_BLACK, 0, 0, 0);
        init_color(NPP_WHITE, 999, 999, 999);
        init_color(NPP_LIGHT_GRAY, 666, 666, 666);
        init_color(NPP_DARK_GRAY, 333, 333, 333);
        init_color(NPP_BROWN, 460, 235, 75);
        init_color(NPP_RED, 999, 0, 0);
        init_color(NPP_MAROON, 450, 0, 0);
        init_color(NPP_ORANGE, 999, 450, 0);
        init_color(NPP_YELLOW, 999, 999, 0);
        init_color(NPP_LIME, 0, 999, 0);
        init_color(NPP_GREEN, 0, 450, 0);
        init_color(NPP_CYAN, 0, 999, 999);
        init_color(NPP_TEAL, 0, 450, 450);
        init_color(NPP_BLUE, 0, 0, 999);
        init_color(NPP_MAGENTA, 999, 0, 999);
        init_color(NPP_PURPLE, 450, 0, 450);
        
        init_pair(0, NPP_WHITE, NPP_BLACK);
        for (unsigned char i = 1; i < 16; i++) {
            init_pair(i, i, NPP_BLACK);
        }

        if (useMouse) {
            printf("\033[?1003h\n");
            mousemask(ALL_MOUSE_EVENTS | REPORT_MOUSE_POSITION, NULL);
        }

        mwin = Window();

        return true;
    }

    /// @brief End - Acts as the ncursespp version of endwin() with some extra cleanup
    /// @param useMouse Whether to clean up the mouse settings or not
    /// @param funcReturn An integer for the function to return (especially useful for reducing clutter in a main() function or smth)
    /// @returns funcReturn
    int end(bool useMouse = false, int funcReturn = 0) {
        endwin();
        // if (useMouse) {printf("\033[?1003h\n");}
        return funcReturn;
    }
}

#endif
