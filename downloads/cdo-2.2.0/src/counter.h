#ifndef COUNTER_H
#define COUNTER_H

namespace cdo
{
class Counter
{
public:
  void start();
  void stop();

  double
  cputime()
  {
    return m_cputime;
  }

private:
  double m_cputime = 0.0;
  char mark[32] = { 0 };
};
}  // namespace cdo

#endif /* COUNTER_H */
