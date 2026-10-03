#ifndef flick_time_point
#define flick_time_point

#include <iostream>
#include <cmath>
#include <iterator>
#include <sstream>
#include <stdexcept>

namespace flick {  
  class time_point
  // https://en.wikipedia.org/wiki/Julian_day
  {
    int year_;
    int month_;
    int day_;
    int hour_;
    int minute_;
    double second_;
  public:
    struct epoch {
      static time_point matlab_datenum() {
	return time_point(-1,12,31,0,0,0);
      }
      static time_point julian_day() {
	return time_point(-4713,11,24,12,0,0);
      }
      static time_point J2000() {
	return time_point(2000,1,1,12,0,0);
      }
    };
    time_point() : time_point(epoch::julian_day()) {
    }
    time_point(int year, int month, int day,
	       int hour, int minute, double second)
      : year_{year}, month_{month}, day_{day},
	hour_{hour}, minute_{minute}, second_{second} {
      ensure(valid());
    }
    double julian_date() const {
      return julian_day_number()+(hour_-12)/24.+minute_/1440.+second_/86400.;
    }
    int year() const {
      return year_;
    }
    int month() const {
      return month_;
    }
    int day() const {
      return day_;
    }
    int hour() const {
      return hour_;
    }
    int minute() const {
      return minute_;
    }
    double second() const {
      return second_;
    }
    double day_of_year() const {
      return julian_date() - time_point(year_,1,1,0,0,0).julian_date();
    }
    double hour_of_day() const {
      double days = julian_date() -
	time_point(year_,month_,day_,0,0,0).julian_date();
      return days*24;
    }
    double J2000() const {
      return julian_date() - epoch::J2000().julian_date();
    }
    double matlab_datenum() const {
      return julian_date() - epoch::matlab_datenum().julian_date();
    }
    friend std::ostream& operator<<(std::ostream& os, const time_point& t) {
      os << t.year_ << " " << t.month_ << " " << t.day_ << " "
	 << t.hour_ << " " << t.minute_ << " " << t.second_;
      return os;
    }
    friend std::istream& operator>>(std::istream& is, time_point& t) {
      int year, month, day, hour, minute;
      double second;
      if (is >> year >> month >> day >> hour >> minute >> second) {
	try {
	  t = time_point{year, month, day, hour, minute, second};
	} catch (const std::exception&) {
	  is.setstate(std::ios::failbit);
	}
      }
      return is;
    }
  private:
    int julian_day_number() const {
      int Y = year_;
      int M = month_;
      int D = day_;
      return (1461 * (Y + 4800 + (M - 14)/12))/4 +
	(367 * (M - 2 - 12 * ((M - 14)/12)))/12 -
	(3 * ((Y + 4900 + (M - 14)/12)/100))/4 + D - 32075;
    }
    bool valid() const {
      return month_ >= 1 && month_ <= 12 &&
	day_ >= 1 && day_ <= 31 &&
	hour_ >= 0 && hour_ <= 23 &&
	minute_ >= 0 && minute_ <= 59 &&
	second_ >= 0 && second_ < 60;
    }
    void ensure(bool b) const {
      if (not b) {
	std::ostringstream os;
	os << "time_point " << *this;
	throw std::runtime_error(os.str());
      }
    }
  };

  inline double datenum_to_julian_date(double matlab_datenum) {
    return matlab_datenum + time_point::epoch::matlab_datenum().julian_date();
  }
  
  inline double J2000_to_julian_date(double J2000) {
    return J2000 + time_point::epoch::J2000().julian_date();
  }

  inline time_point make_time_point(double julian_date)
  // https://en.wikipedia.org/wiki/Julian_day
  {
    int y = 4716;
    int j = 1401;
    int m = 2;
    int n = 12;
    int r = 4;
    int p = 1461;
    int v = 3;
    int u = 5;
    int s = 153;
    int w = 2;
    int B = 274277;
    int C = -38;      
    int J = round(julian_date);
    int f = J+j+(((4*J+B)/146097)*3)/4+C;
    int e = r*f+v;
    int g = (e % p)/r;
    int h = u*g+w;
    int D = (h % s)/u+1;
    int M = (h/s+m) % n + 1;
    int Y = e/p-y+(n+m-M)/n;
    int hr = (julian_date-time_point(Y,M,D,0,0,0).julian_date())*24;
    int mi = (julian_date-time_point(Y,M,D,hr,0,0).julian_date())*60*24;
    double se = (julian_date-time_point(Y,M,D,hr,mi,0).julian_date())*60*60*24;
    return time_point(Y,M,D,hr,mi,se);
  }
  class time_converter {
    struct leap_second {
      int year;
      int month;
      int day;
      int tai_minus_utc;
    };

    static constexpr leap_second leap_seconds_[] = {
      {1972, 1, 1, 10}, {1972, 7, 1, 11}, {1973, 1, 1, 12},
      {1974, 1, 1, 13}, {1975, 1, 1, 14}, {1976, 1, 1, 15},
      {1977, 1, 1, 16}, {1978, 1, 1, 17}, {1979, 1, 1, 18},
      {1980, 1, 1, 19}, {1981, 7, 1, 20}, {1982, 7, 1, 21},
      {1983, 7, 1, 22}, {1985, 7, 1, 23}, {1988, 1, 1, 24},
      {1990, 1, 1, 25}, {1991, 1, 1, 26}, {1992, 7, 1, 27},
      {1993, 7, 1, 28}, {1994, 7, 1, 29}, {1996, 1, 1, 30},
      {1997, 7, 1, 31}, {1999, 1, 1, 32}, {2006, 1, 1, 33},
      {2009, 1, 1, 34}, {2012, 7, 1, 35}, {2015, 7, 1, 36},
      {2017, 1, 1, 37}
    };

    static double TT_minus_UTC(const time_point& utc) {
      double jd = utc.julian_date();
      for (auto it = std::rbegin(leap_seconds_);
           it != std::rend(leap_seconds_); ++it) {
        if (jd >= time_point(it->year, it->month, it->day, 0, 0, 0).julian_date())
          return it->tai_minus_utc + 32.184;
      }
      throw std::domain_error("UTC to TT conversion is supported from 1972-01-01");
    }

    static double TT_minus_UTC(double tt_julian_date) {
      for (auto it = std::rbegin(leap_seconds_);
           it != std::rend(leap_seconds_); ++it) {
        double utc_boundary =
          time_point(it->year, it->month, it->day, 0, 0, 0).julian_date();
        double tt_boundary = utc_boundary + (it->tai_minus_utc + 32.184)/86400;
        if (tt_julian_date >= tt_boundary)
          return it->tai_minus_utc + 32.184;
      }
      throw std::domain_error("TT to UTC conversion is supported from 1972-01-01");
    }
  public:
    time_point UTC_to_TT(const time_point& utc) {
      return make_time_point(utc.julian_date()+TT_minus_UTC(utc)/86400);
    }
    time_point TT_to_UTC(const time_point& tt) {
      return make_time_point(tt.julian_date()-TT_minus_UTC(tt.julian_date())/86400);
    }
  };
}

#endif
