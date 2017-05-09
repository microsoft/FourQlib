#ifndef STM32F4_WRAPPER_H
#define STM32F4_WRAPPER_H


void clock_setup(void);
void gpio_setup(void);
void usart_setup(int baud);
void rng_setup(void);
void signal_host(void);
void send_USART_str(const unsigned char* in);
void random_int(uint32_t*,int);


#endif
