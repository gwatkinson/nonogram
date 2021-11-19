import numpy as np
import matplotlib.pyplot as plt
from IPython.display import display, clear_output

def generate_spaces(p, nc):
# def generate_spaces(number_of_spaces_int, number_of_coef_int):
    """
    Generate all possible configuration of `n` spaces between `n` coefficients. 
    
    Parameters
    ----------
    number_of_space_int : int
        Number of spaces in total to distribute between the coefficients.
    number_of_coef_int : int
        Number of coefficients.
    
    Note
    ----
    This function is recursive. And the number of configuration is the Pascal coefficients of
    `number_of_coef_int` in `number_of coef_int + number_of_spaces_int`.
    
    Returns
    -------
    list
        List of configurations of given spaces between the given number of coefficients.
        
        Each configuration is a list of size `number_of_coef_int+1` and the sum of 
        the elements in the list is equal to `number_of_spaces_int`.
    """
    # options_list = []
    # 
    # if number_of_coef_int == 0:
    #     return [[number_of_spaces_int]]
    # else:
    #     for first_space_int in range(number_of_spaces_int + 1):
    #         number_of_coef_int_recursive = number_of_coef_int - 1
    #         number_of_spaces_int_recursive = number_of_spaces_int - first_space_int
    #         options_list_recursive = generate_spaces(number_of_spaces_int_recursive, number_of_coef_int_recursive)
    #         for option_list in options_list_recursive:
    #             tmp_list = [first_space_int]
    #             tmp_list.extend(option_list)
    #             options_list.append(tmp_list)
    #     return options_list
    options = []
    
    if nc == 0:
        return [[p]]
    else:
        for x in range(p+1):
            nc_p = nc - 1
            p_p = p - x
            options_p = generate_spaces(p_p, nc_p)
            for o in options_p:
                l = [x]
                l.extend(o)
                options.append(l)
        return options 
    
    
def generate_options(n, coefs):
# def generate_options(length_of_row_int, coefficients_list):
    """
    Generate all possible configuration for a given row with the coefficients. 
    
    Parameters
    ----------
    length_of_row_int : int
        Length of the row to fill.
    number_of_coef_int : int
        List of the coefficients for the row.
    
    Note
    ----
    This function uses :func:`generate_space`. It is used to initialise the solver.
    
    Returns
    -------
    list
        List of all possible options for the repartion of the coefficients for the specified row.
        
        An option is a list of `0` and `1`, with `1` meaning the presence of a square and 
        `0` the presence of a cross. Every options is of the length of `length_of_row_int` and
        is valid (i.e. coefficients are all here).
    
    Raises
    ------
    ValueError
        If there is no configuration possible, when the sum of the coefficients is to big for the row.
    """
    # number_of_coef_int = len(coefficients_list)
    # number_of_spaces_int = length_of_row_int - sum(coefficients_list) + number_of_coef_int - 1
    # 
    # if number_of_spaces_int < 0:
    #     raise ValueError(f"The sum of the coefficients is too big. \n"
    #     + "Should be less than {lenght_of_row_int}, actually is {lenght_of_row_int - number_of_spaces_int}")
    # 
    # options_list = []
    # space_configurations_list = generate_spaces(number_of_spaces_int, number_of_coef_int)
    # 
    # for space_list in space_configurations_list:
    #     row_list = [0 for i in range(length_of_row_int)]
    #     position_int = 0
    #     for i in range(len(space_list) - 1):
    #         logging.info("Start")
    #         position_int += space_list[i]
    #         logging.info(f"Space: {space_list[i]}")
    #         logging.info(f"Position: {position_int}")
    #         logging.info(f"Coef:{coefficients_list[i]}")
    #         logging.info(f"Position: {position_int + coefficients_list[i]}")
    #         for j in range(coefficients_list[i]):
    #             row_list[position_int + j] = 1
    #                 
    #         position_int += coefficients_list[i] + 1
    #     options_list.append(row_list)
    # return options_lis
    nc = len(coefs)
    t = sum(coefs) + nc - 1
    p = n-t
    
    if t > n:
        raise ValueError(f"The sum of the coefficients is too big. \n"
        + "Should be less than {n}, actually is {t}")
    
    options = []
    spaces = generate_spaces(p, nc)
    
    for s in spaces:
        l = [0 for i in range(n)]
        res = 0
        for i in range(len(s)-1):
            res += s[i]
            for j in range(coefs[i]):
                l[res+j] = 1
            res += coefs[i] + 1
        options.append(l)
    return options

class Game:
    def __init__(self, board, row_coefs, col_coefs):
        """
        Game class containing the board, the coefficients of the rows and columns and a display function.
        
        Parameters
        ----------
        board : Array
            Numpy array containing the initial game state.
        row_coefs : list[list[int]]
            Contains the list of cofficients for the rows.
        col_coefs : list[list[int]]
            Contains the list of cofficients for the columns.
        
        Attributes
        ----------
        board : Array
            The board of the game.
        n, row_length : int
            The length of the rows.
        m, col_length : int
            The length of the columns.
        row_coefs : list
            The coefficients for the rows.
        col_coefs : list
            The coefficients for the columns.
        fig : Figure
            Matplotlib figure.
        ax : Axis
            Matplotlib axis.
        """
        self.board = board
        self.row_length = self.n = board.shape[0]
        self.col_length = self.m = board.shape[1]
        self.row_coefs = row_coefs
        self.col_coefs = col_coefs
        self.fig, self.ax = plt.subplots()
        self.assert_game()
    
    def display(self, iteration=None):
        """Dipslay the board as a grayscale heatmap."""
        self.ax.clear()
        self.ax.imshow(self.board, cmap='gray_r', vmin=-1, vmax=1)
        # self.ax.tick_params(
        #     axis='both',
        #     which='both',
        #     bottom='off',
        #     top='off',
        #     labelbottom='off',
        #     right='off',
        #     left='off',
        #     labelleft='off'
        # )
        if iteration:
            self.ax.set_title(f"Iteration: {iteration}")
        
    def assert_game(self):
        """Assert that the given initial values are valid."""
        assert_rows_list = [i for i in range(len(self.row_coefs)) if (sum(self.row_coefs[i]) + len(self.row_coefs[i]) - 1 > self.n)]
        assert_cols_list = [i for i in range(len(self.col_coefs)) if (sum(self.col_coefs[i]) + len(self.col_coefs[i]) - 1 > self.m)]
        assert_values_list = np.argwhere(~np.isin(self.board, [-1, 0, 1])).tolist()
        
        assert_errors = {
            "rows": (assert_rows_list, f"Coefficients for rows {assert_rows_list} are invalid"),
            "cols": (assert_cols_list, f"Coefficients fo cols {assert_cols_list} are invalid"),
            "board_values": (assert_values_list, f"Values {assert_values_list} are invalid, should be in `[-1, 0, 1]`")
        }
        assert_message = "\n".join([assert_errors[key][1] for key in assert_errors if assert_errors[key][0]])
        
        if assert_message:
            raise AssertionError(assert_message)
            
            
class Solver:
    def __init__(self, game):
        """
        Solver class that contains methods to solve the game.
        
        Parameters
        ----------
        game : :obj:`Game`
            The ``Game`` class object to solve.
        
        Attributes
        ----------
        game : :obj:`Game`
            The ``Game`` class.game.
        board : Array
            The board of the game.
        n, row_length : int
            The length of the rows.
        m, col_length : int
            The length of the columns.
        row_coefs : list
            The coefficients for the rows.
        col_coefs : list
            The coefficients for the columns.
        row_options : list
            The list of the options for the rows.    
            It is a list of list of list of int.    
            It is initialised with the ``init_options("row")`` function.
        col_options : list
            The list of the options for the columns.    
            It is a list of list of list of int.    
            It is initialised with the ``init_options("col")`` function.
        """
        self.game = game
        self.board = game.board.copy()
        self.previous_board = None
        self.row_length = self.n = game.n
        self.col_length = self.m = game.m
        self.row_coefs = game.row_coefs
        self.col_coefs = game.col_coefs
        self.row_options = self.init_options("row")
        self.col_options = self.init_options("col")
    
    
    def init_options(self, axis):
        """Generate all the option for the given axis with ``generate_options``."""
        if axis == "row" or axis =="r":
            size, list_coefs = self.row_length, self.row_coefs
        elif axis == "col" or axis =="c":
            size, list_coefs = self.col_length, self.col_coefs
        else:
            raise ValueError("Parameter must be in 'row', 'r' for rows or 'col', 'c' for columns")
        
        options = []
        for i in range(len(list_coefs)):
            options.append(generate_options(size, list_coefs[i]))
        
        return options
    
    
    def is_finished(self):
        """Checks if the game is finished, by checking that there are no `-1` left."""
        return (self.board != -1).all()
    
    
    def check_absurd_options(self):
        """Remove options that are absurd given the current board."""
        # Rows
        for r in range(self.n):
            row = self.board[r, :]
            for x, opt in enumerate(self.row_options[r]):
                if np.any([(row[i]!=-1 and opt[i]!=row[i]) for i in range(len(row))]):
                    self.row_options[r].pop(x)
        
        # Columns
        for c in range(self.m):
            col = self.board[:, c]
            for x, opt in enumerate(self.col_options[c]):
                if np.any([(col[i]!=-1 and opt[i]!=col[i]) for i in range(len(col))]):
                    self.col_options[c].pop(x)
                    
        return None
    
    
    def check_superposition(self):
        """
        Fill the board with the superposition method.
        
        Description of the method on `nonogram.org`_.
        
        .. _nonogram.org:
            https://www.nonograms.org/methods
        """
        # Check rows
        for r in range(self.n):
            for i in range(self.m):
                if self.board[r, i] == -1 and (np.array(self.row_options[r])[:, i] == self.row_options[r][0][i]).all():
                    self.board[r, i] = self.row_options[r][0][i]
        
        # Check columns
        for c in range(self.m):
            for i in range(self.n):
                if self.board[i, c] == -1 and (np.array(self.col_options[c])[:, i] == self.col_options[c][0][i]).all():
                    self.board[i, c] = self.col_options[c][0][i]
                    
        return None
    
    
    def display_board(self, fig, ax, iteration=None):
        """Dipslay the board as a grayscale heatmap."""
        ax.cla()
        ax.imshow(self.board, cmap='gray_r', vmin=-1, vmax=1)
        if iteration:
            ax.set_title(f"Iteration: {iteration}")
        display(fig)
        clear_output(wait = True)
          
    
    def solve(self, check_previous=True, display=False, pause=0.0):
        """
        Solve the game by iterating and checking for absurd values and superposition to fill the board.
        
        Parameters
        ----------
        check_previous : bool
            True if you check that there is progress every iteration.    
            Stops it if there is no more progress.
        display : bool
            True for displaying progress every iterations.
        pause : float
            Time to pause after displaying an iteration.
        """
        if display:
            fig, ax = plt.subplots()
        i = 0
        while not self.is_finished():
            if display:
                self.display_board(fig, ax, iteration=i)
                clear_output(wait = True)
                if pause > 0:
                    plt.pause(0.5)
            if check_previous:
                self.previous_board = self.board.copy()
            
            self.check_absurd_options()
            self.check_superposition()
            
            i += 1
            
            if check_previous and (self.previous_board == self.board).all():
                raise Error(f"No Progress. Stopped at iteration {i}")
            
        if self.is_finished():
            print(f"Game is solved in {i} iterations!\n")
            self.display_board(fig, ax, iteration=i)