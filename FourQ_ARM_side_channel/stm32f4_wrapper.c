#include <stdint.h>
#include <rcc.h>
#include <gpio.h>
#include <usart.h>
#include <nvic.h>
#include <rng.h>
#include "stm32f4_wrapper.h"


void clock_setup(void)
{
	rcc_clock_setup_hse_3v3(&hse_8mhz_3v3[CLOCK_3V3_48MHZ]);
	rcc_periph_clock_enable(RCC_GPIOD);
    rcc_periph_clock_enable(RCC_GPIOA);
    rcc_periph_clock_enable(RCC_USART2);
    rcc_periph_clock_enable(RCC_DMA1);
	rcc_periph_clock_enable(RCC_RNG);
}

void gpio_setup(void)
{
    gpio_mode_setup(GPIOA, GPIO_MODE_AF, GPIO_PUPD_NONE, GPIO2 | GPIO3);
    gpio_set_af(GPIOA, GPIO_AF7, GPIO2 | GPIO3);
	gpio_mode_setup(GPIOD, GPIO_MODE_OUTPUT,GPIO_PUPD_NONE, GPIO12);
}

void usart_setup(int baud)
{
    usart_set_baudrate(USART2, baud);
    usart_set_databits(USART2, 8);
    usart_set_stopbits(USART2, USART_STOPBITS_1);
    usart_set_mode(USART2, USART_MODE_TX_RX);
    usart_set_parity(USART2, USART_PARITY_NONE);
    usart_set_flow_control(USART2, USART_FLOWCONTROL_NONE);
    usart_enable(USART2);
	gpio_set(GPIOD, GPIO12);
}

void rng_setup(void)
{
	RNG_CR |= RNG_CR_IE;
	RNG_CR |= RNG_CR_RNGEN;
}

void send_USART_str(const unsigned char* in)
{
    gpio_toggle(GPIOD, GPIO12);
	int i;

    for(i = 0; in[i] != 0; i++) {
		usart_send_blocking(USART2, in[i]);
    }    
	usart_send_blocking(USART2, '\r');
    usart_send_blocking(USART2, '\n');
}

void signal_host(void) 
{
    usart_send_blocking(USART2, (char)4);
}

void random_int(uint32_t* urnd, int n)
{
	unsigned int last_value=0;
	unsigned int new_value=0;
	int i;
	unsigned int error_bits = 0;

	for(i = 0; i < n; i++)
	{
	    error_bits = RNG_SR_SEIS | RNG_SR_CEIS;
	    while (new_value == last_value) {
		    if (((RNG_SR & error_bits) == 0) && ((RNG_SR & RNG_SR_DRDY) == 1)) {
			    new_value = RNG_DR;
		    }
	    }
	    last_value = new_value;
	    urnd[i] = new_value;
	}
}