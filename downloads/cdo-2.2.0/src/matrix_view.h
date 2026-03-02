#ifndef MATRIX_VIEW_H
#define MATRIX_VIEW_H

// Modified code from https://github.com/pwwiur/Matrix

template <typename T>
class MatrixView
{
private:
  T *arr;
  size_t numRows, numColumns;

  class Proxy
  {
  public:
    Proxy(T *arr, size_t columns, size_t x) : px_arr(arr), px_columns(columns), px_x(x) {}
    T &
    operator[](size_t y)
    {
      return px_arr[y + px_columns * px_x];
    }

  private:
    T *px_arr;
    size_t px_columns, px_x;
  };

public:
  MatrixView(T *array, size_t rows, size_t columns) : arr(array), numRows(rows), numColumns(columns)
  {
    if (rows == 0 || columns == 0) throw "Matrix initialization arguments rows and columns must be greater than zero.";
  }

  Proxy
  operator[](size_t x)
  {
    return Proxy(arr, numColumns, x);
  }
};

#endif  // MATRIX_VIEW_H
