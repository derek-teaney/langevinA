#ifndef MODELA_MAKE_UNIQUE_H
#define MODELA_MAKE_UNIQUE_H


template<typename T, typename... Args>
std::unique_ptr<T> make_unique(Args&&... args)
{
#ifndef MODELA_NO_MAKE_UNIQUE
  return std::make_unique<T>(std::forward<Args>(args)...);
#else 
  return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
#endif
}
#endif


