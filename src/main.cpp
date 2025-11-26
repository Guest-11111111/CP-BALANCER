#include "vkey/vkey.hpp"
#include <appdef.h>
#include <sdk/os/debug.h>
#include <sdk/os/input.h>
#include <sdk/os/lcd.h>
#include <sdk/os/string.h>

#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <cctype>

APP_NAME("CP-Balancer")
APP_AUTHOR("Guest-11111111")
APP_DESCRIPTION("Classpad-HHK3 Chemical Equation Balancer")
APP_VERSION("2.0.0")

#define MAX_EQU_LEN 128

// --- Utility Functions ---
unsigned int str_len(const char* s) { return String_Strlen(s); }
char* str_cat(char* d, const char* s) { return String_Strcat(d, s); }
char* str_cpy(char* d, const char* s) { return String_Strcpy(d, s); }
int str_cmp(const char* a, const char* b) { return String_Strcmp(a, b); }
const char* str_chr(const char* s, char c) { return String_Strchr(s, c); }

void trim_spaces(char* s) {
    int n = str_len(s);
    int j = 0;
    for (int i = 0; i < n; ++i) {
        if (!isspace(s[i])) s[j++] = s[i];
    }
    s[j] = 0;
}

// --- Input Functions ---
void get_input_equation(char* buf, int buflen) {
    int pos = 0;
    Debug_SetCursorPosition(0, 0);
    Debug_PrintString("Enter equation (EXE to finish):", false);
    LCD_Refresh();
    while (pos < buflen - 1) {
        struct Input_Event ev = {};
        GetInput(&ev, 0xFFFFFFFF, 0x10);
        if (ev.type == EVENT_KEY && ev.data.key.direction == KEY_PRESSED) {
            int kc = ev.data.key.keyCode;
            if (kc == KEYCODE_EXE) break;
            if (kc == KEYCODE_BACKSPACE && pos > 0) { pos--; }
            else if (kc == KEYCODE_7) buf[pos++] = 'H';
            else if (kc == KEYCODE_8) buf[pos++] = 'I';
            else if (kc == KEYCODE_9) buf[pos++] = 'J';
            else if (kc == KEYCODE_PLUS) buf[pos++] = '+';
            else if (kc == KEYCODE_EQUALS) buf[pos++] = 'A';
            else if (kc == KEYCODE_MINUS) buf[pos++] = '-';
            else if (kc == KEYCODE_COMMA) buf[pos++] = ',';
            else if (kc == KEYCODE_OPEN_PARENTHESIS) buf[pos++] = 'G';
            else if (kc == KEYCODE_CLOSE_PARENTHESIS) buf[pos++] = ')';
            else if (kc == KEYCODE_X) buf[pos++] = 'B';
            else if (kc == KEYCODE_Y) buf[pos++] = 'C';
            else if (kc == KEYCODE_Z) buf[pos++] = 'D';
            else if (kc == KEYCODE_DIVIDE) buf[pos++] = 'F';
            else if (kc == KEYCODE_POWER) buf[pos++] = 'E';
            else if (kc == KEYCODE_TIMES) buf[pos++] = 'K';
            else if (kc == KEYCODE_SHIFT) buf[pos++] = 'L';
        }
        buf[pos] = 0;
        Debug_SetCursorPosition(1, 0);
        Debug_PrintString(buf, false);
        LCD_Refresh();
    }
    buf[pos] = 0;
}

// --- Parse formula to element counts ---
void parse_formula(const char* formula, std::map<std::string, int>& elem_counts) {
    int i = 0, n = str_len(formula);
    while (i < n) {
        if (isupper(formula[i])) {
            char elem[3] = { formula[i], 0, 0 };
            if (i+1 < n && islower(formula[i+1])) {
                elem[1] = formula[i+1];
                ++i;
            }
            int count = 0;
            int j = i+1;
            while (j < n && isdigit(formula[j])) {
                count = count * 10 + (formula[j] - '0');
                ++j;
            }
            if (count == 0) count = 1;
            elem_counts[std::string(elem)] += count;
            i = j-1;
        }
        ++i;
    }
}

// --- Split equation into lhs/rhs formulas ---
bool split_equation(const char* eq, std::vector<std::string>& lhs, std::vector<std::string>& rhs) {
    const char* eq_pos = str_chr(eq, '=');
    if (!eq_pos) return false;
    int left_len = eq_pos - eq;
    char left[MAX_EQU_LEN] = {0}, right[MAX_EQU_LEN] = {0};
    for (int i = 0; i < left_len; ++i) left[i] = eq[i];
    left[left_len] = 0;
    str_cpy(right, eq_pos + 1);
    trim_spaces(left);
    trim_spaces(right);

    auto split_side = [](const char* str, std::vector<std::string>& out) {
        int n = str_len(str);
        int last = 0;
        for (int i = 0; i <= n; ++i) {
            if (str[i] == '+' || str[i] == 0) {
                char buf[32] = {0};
                int len = i - last;
                if (len > 0 && len < 31) {
                    for (int j = 0; j < len; ++j) buf[j] = str[last + j];
                    buf[len] = 0;
                    out.push_back(std::string(buf));
                }
                last = i + 1;
            }
        }
    };
    split_side(left, lhs);
    split_side(right, rhs);
    return true;
}

// --- Build matrix for balancing ---
void build_matrix(const std::vector<std::string>& lhs, const std::vector<std::string>& rhs,
                  std::vector<std::string>& elements,
                  std::vector<std::vector<int> >& matrix) {
    std::map<std::string, int> elem_map;
    std::vector<std::map<std::string, int> > all_counts;

    for (size_t i = 0; i < lhs.size(); ++i) {
        std::map<std::string, int> counts;
        parse_formula(lhs[i].c_str(), counts);
        all_counts.push_back(counts);
        for (auto& p : counts) elem_map[p.first] = 1;
    }
    for (size_t i = 0; i < rhs.size(); ++i) {
        std::map<std::string, int> counts;
        parse_formula(rhs[i].c_str(), counts);
        all_counts.push_back(counts);
        for (auto& p : counts) elem_map[p.first] = 1;
    }
    elements.clear();
    for (auto& p : elem_map) elements.push_back(p.first);

    int rows = elements.size();
    int cols = lhs.size() + rhs.size();
    matrix.assign(rows, std::vector<int>(cols, 0));
    for (int i = 0; i < rows; ++i) {
        for (size_t j = 0; j < lhs.size(); ++j)
            matrix[i][j] = all_counts[j][elements[i]];
        for (size_t j = 0; j < rhs.size(); ++j)
            matrix[i][lhs.size() + j] = -all_counts[lhs.size() + j][elements[i]];
    }
}

// --- GCD for normalization ---
int gcd(int a, int b) {
    while (b) { int t = b; b = a % b; a = t; }
    return a < 0 ? -a : a;
}
int vector_gcd(const std::vector<int>& v) {
    int res = 0;
    for (size_t i = 0; i < v.size(); ++i) if (v[i]) res = res ? gcd(res, v[i]) : v[i];
    return res ? res : 1;
}

// --- Gaussian elimination for integer solution ---
bool solve_matrix(const std::vector<std::vector<int> >& matrix, std::vector<int>& coeffs) {
    int rows = matrix.size(), cols = matrix[0].size();
    std::vector<std::vector<int>> mat = matrix;
    std::vector<int> solution(cols, 0);
    // Augment with zero column for homogeneous system
    for (int i = 0; i < rows; ++i) mat[i].push_back(0);
    // Gaussian elimination
    int rank = 0;
    for (int col = 0; col < cols; ++col) {
        int pivot = -1;
        for (int row = rank; row < rows; ++row) {
            if (mat[row][col] != 0) { pivot = row; break; }
        }
        if (pivot == -1) continue;
        std::swap(mat[rank], mat[pivot]);
        int div = mat[rank][col];
        for (int j = col; j <= cols; ++j) mat[rank][j] /= gcd(div, mat[rank][j]);
        for (int row = 0; row < rows; ++row) {
            if (row != rank && mat[row][col] != 0) {
                int factor = mat[row][col];
                for (int j = col; j <= cols; ++j)
                    mat[row][j] -= factor * mat[rank][j] / mat[rank][col];
            }
        }
        ++rank;
    }
    // Find a nontrivial integer solution
    solution[cols-1] = 1;
    for (int i = rows-1; i >= 0; --i) {
        int s = 0, idx = -1;
        for (int j = 0; j < cols; ++j) {
            if (mat[i][j] != 0) {
                if (idx == -1) idx = j;
                else s += mat[i][j] * solution[j];
            }
        }
        if (idx != -1) solution[idx] = -s / mat[i][idx];
    }
    int d = vector_gcd(solution);
    for (size_t i = 0; i < solution.size(); ++i) solution[i] /= d;
    coeffs = solution;
    // Check for all positive
    bool allpos = std::all_of(coeffs.begin(), coeffs.end(), [](int x){return x>0;});
    if (!allpos) {
        int sign = 1;
        for (int x : coeffs) if (x < 0) sign = -1;
        for (size_t i = 0; i < coeffs.size(); ++i) coeffs[i] *= sign;
    }
    return true;
}

void display_balanced(const std::vector<std::string>& lhs, const std::vector<std::string>& rhs, const std::vector<int>& coeffs) {
    char buf[MAX_EQU_LEN] = {0};
    for (size_t i = 0; i < lhs.size(); ++i) {
        if (coeffs[i] != 1) {
            char tmp[8]; snprintf(tmp, 8, "%d", coeffs[i]); str_cat(buf, tmp);
        }
        str_cat(buf, lhs[i].c_str());
        if (i < lhs.size()-1) str_cat(buf, " + ");
    }
    str_cat(buf, " = ");
    for (size_t i = 0; i < rhs.size(); ++i) {
        if (coeffs[lhs.size()+i] != 1) {
            char tmp[8]; snprintf(tmp, 8, "%d", coeffs[lhs.size()+i]); str_cat(buf, tmp);
        }
        str_cat(buf, rhs[i].c_str());
        if (i < rhs.size()-1) str_cat(buf, " + ");
    }
    Debug_Printf(1,1,false, 0, "Balanced: %s", buf);
    LCD_Refresh();
    Debug_WaitKey();
}

void show_error(const char* msg) {
    Debug_SetCursorPosition(0, 0);
    Debug_PrintString(msg, false);
    LCD_Refresh();
    Debug_WaitKey();
}

int main() {
    LCD_ClearScreen();
    LCD_Refresh();
    char equation[MAX_EQU_LEN] = {0};
    get_input_equation(equation, MAX_EQU_LEN);

    std::vector<std::string> lhs, rhs;
    if (!split_equation(equation, lhs, rhs)) {
        show_error("Invalid format! Use e.g. H2+O2=H2O");
        return 1;
    }
    std::vector<std::string> elements;
    std::vector<std::vector<int> > matrix;
    build_matrix(lhs, rhs, elements, matrix);

    std::vector<int> coeffs;
    if (!solve_matrix(matrix, coeffs)) {
        show_error("Could not balance equation.");
        return 1;
    }
    display_balanced(lhs, rhs, coeffs);
    return 0;
}
