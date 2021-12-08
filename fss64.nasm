%include "sseutils64.nasm"
; min_vector_64(VECTOR x, int n, type* max)

section .data
section .bss
section .text

global min_vector_64
	; rdi = x
    ; rsi = n
    ; rdx = max
    UNROLL_MAX equ 8 
min_vector_64 :
	start

	mov     rax,rdi                 ; x
	mov     rdi,rsi                 ; n
	sub     rdi,UNROLL_MAX-1        ; gestione vettore non multiplo
	vmovapd ymm0,[rax]              ; primi quattro elementi (in teoria)
	mov 	rsi, 4                  ; i=4

fori_min:
    vmovapd ymm1,[rax+rsi*8]        ; xmm1<-x[...]
	vminpd 	ymm0,ymm1               ; confronto ymm0 ymm1

    vmovapd ymm1,[rax+rsi*8+32]     ; UNROLL    
	vminpd 	ymm0,ymm1

	add	    rsi,UNROLL_MAX          ; i+=UNROLL

	cmp 	rsi,rdi                 ; i<n-7?
	jl 	    fori_min


	add 	rdi,UNROLL_MAX-1        ; ripristino n



                                            ; riduzione vettore
	vmovapd     ymm1,ymm0                   ; copia su ymm1
	vperm2f128  ymm0,ymm0,ymm0,00110000b    ; permutazione dei 128 bit
	vminpd 	    ymm0,ymm1                   
    vmovapd     ymm1,ymm0                   ; copia su ymm1
    vshufpd     xmm1,xmm1,01b               ; per confrontare gli ultimi due double si passa
	vminpd 	    xmm0,xmm1                   ; agli xmm e si usa shufpd per invertire i 64 bit

	cmp 	rsi,rdi                 ; i<n? 
	jge 	end_min                 ; gestione caso lunghezza multipla di 8

forino_min:                         ; replica codice per versione scalare double
	vmovq 	xmm1,[rax+rsi*8]    	; muovo un double in xmm1. Nota : la v sembrerebbe non necessaria ma dovrebbe impedire il context switch evitando di incorrere nella penalità    
	vminpd 	xmm0, xmm1				; la v è presente per lo stesso ragionamento di sopra
	add	    rsi,1
    cmp     rsi,rdi
	jl	    forino_min

end_min:                            ; caricamento del minimo in memoria
	vmovq 	[rdx], xmm0			

	stop
	

;vector_sum_64(MATRIX x, int offset, int n,VECTOR v);

section .data
section .bss
section .text

global vector_sum_64
    ; rdi = x
    ; rsi = offset
    ; rdx = n
    ; rcx = v

    UNROLL_VS equ 8 

	; controllare sopra 

vector_sum_64:
    start

    mov     rax,rdi                 ; x
    mov     rbx,rsi                 ; offset
    imul    rbx,8                   ; porta offset a versione byte. Un double è 8 byte
    add     rax,rbx                 ; porta l'indice alla posizione del vettore target

    mov     rbx,rcx                 ; v

    mov     rdi,rdx                 ; n
    sub     rdi,UNROLL_VS-1         ; unroll

    mov     rsi,0                   ; i=0
fori_vs:

    vmovapd  ymm0,[rbx+rsi*8]        ; v[...]
    vaddpd   ymm0,[rax+rsi*8]        ; somma v[...] con x[...]
    vmovapd  [rax+rsi*8],ymm0        ; carica il risultato su x[...] in memoria

    vmovapd  ymm0,[rbx+rsi*8+32]     ; UNROLL
    vaddpd   ymm0,[rax+rsi*8+32]
    vmovapd  [rax+rsi*8+32],ymm0

    add     rsi,UNROLL_VS           ; i+=8

    cmp     rsi,rdi                 ; i<n-7?
    jl      fori_vs


    add     rdi,UNROLL_VS-1         ; ripristino n

    cmp     rsi,rdi                 ; i<n?
    jge     end_vs
    
forino_vs:                          ; gestione caso vettore non multiplo di 8
    vmovq   xmm0,[rbx+rsi*8]  
    vmovq   xmm1, [rax+rsi*8] 		; vedere se si possa evitare lo spostamento nel registro e farlo direttamente dalla memoria  
    vaddpd  xmm0,xmm1 
    vmovq   [rcx+rsi*8],xmm0

    inc     rsi
    cmp     rsi,rdi
    jl      forino_vs

end_vs:
    stop


;euclidian_distance_64(MATRIX x, int offset, VECTOR y, int d,type* dist)
section .data
section .bss
section .text

global euclidian_distance_64
    ; rdi = x
    ; rsi = offset
    ; rdx = y
    ; rcx = d
    ; r8 = dist

    UNROLL_EUC equ 8

euclidian_distance_64:
    start
    
    mov     rax,rdi             ; x
    mov     rbx,rsi             ; offset
    mov     rdi,rcx             ; d
    mov     rcx,rdx             ; y
    imul    rbx,8               ; offset in versione byte
    add     rax,rbx             ; indirizzo del vettore target nella matrice
    
    sub     rdi,UNROLL_EUC-1    ; gestione vettore non multiplo di 8
    
    vxorpd  xmm0,xmm0,xmm0          ; ret=0
    mov     rsi,0               ; i=0
fori_euc:   
    vmovapd ymm1,[rax+rsi*8]     ; xmm1<- x[...]
    vmovapd ymm2,[rcx+rsi*8]     ; xmm2<- x[...]
    vsubpd  ymm1,ymm2            ; (x-y)
    vmulpd  ymm1,ymm1            ; (x-y)^2
    vaddpd  ymm0,ymm1            ; ret+=(x-y)^2

    vmovapd ymm1,[rax+rsi*8+32]  ; UNROLL
    vmovapd xmm2,[rcx+rsi*8+32]
    vsubpd  ymm1,ymm2
    vmulpd  ymm1,ymm1
    vaddpd  ymm0,ymm1
    
    add    rsi,UNROLL_EUC       ; i+=8
    cmp    rsi,rdi              ; i<n-7?
    jl     fori_euc
    
    add    rdi,UNROLL_EUC-1        ; i+=8
    
    vhaddpd     ymm0,ymm0                   ; riduzione di ymm0
    vperm2f128  ymm2,ymm0,ymm0,00000011b
    vaddsd      xmm0,xmm2
        
    cmp   rsi,rdi               ; caso vettore multiplo di 8
    jge   end_euc
    
forino_euc:                     ; replica scalare della sezione precedente
    vmovq xmm1,[rax+rsi*8]
    vmovq xmm2,[rcx+rsi*8]
    vsubpd xmm1,xmm2
    vmulpd xmm1,xmm1
    vaddpd xmm0,xmm1
    add rsi,1
    cmp rsi,rdi
    jl forino_euc
    
end_euc:                        
    vsqrtpd xmm0,xmm0            ; sqrt(ret)
    
    vmovq   [r8],xmm0
    
    stop


;eval_f_64(VECTOR x,int d,VECTIOR c, int offset,type* quad,type* scalar)

section .data
section .bss
section .text
global main
    ; rdi = x
    ; rsi = d
    ; rdx = c
    ; rcx = offset
    ; r8 = quad
    ; r9 = scalar
    UNROLL_EVALF equ 8

main:

    start

        mov     rax,rdi
        imul    rcx,8
        add     rax,rcx     ; x[i]

        mov     rbx,rdx     ; c
        mov     rdi,rsi     ; d

        vxorps  ymm0,ymm0   ; quad=0    
        vxorps  ymm1,ymm1   ; scalar=0

        sub     rdi,UNROLL_EVALF-1  ; d-7

        mov     rsi,0               ; i=0
        
fori_evalf:
    vmovapd ymm2,[rax+rsi*8]        
    vmovapd ymm3,ymm2
    vmulpd  ymm3,ymm3
    vaddpd  ymm0,ymm3
    

    vmovapd ymm3,[rbx+rsi*8]
    vmulpd  ymm2,ymm3
    vaddpd  ymm1,ymm2


    vmovapd ymm2,[rax+rsi*8+32]
    vmovapd ymm3,ymm2
    vmulpd  ymm3,ymm3
    vaddpd  ymm0,ymm3
    

    vmovapd ymm3,[rbx+rsi*8+32]
    vmulpd  ymm2,ymm3
    vaddpd  ymm1,ymm2

    add     rsi,UNROLL_EVALF 
    cmp     rsi,rdi
    jl      fori_evalf


    add         rdi,UNROLL_EVALF-1

    vhaddpd     ymm0,ymm0
    vperm2f128  ymm2,ymm0,ymm0,00000011b
    vaddsd      xmm0,xmm2
    


    vhaddpd     ymm1,ymm1
    vperm2f128  ymm2,ymm1,ymm1,00000011b
    vaddsd      xmm1,xmm2

    cmp         rsi,rdi
    jge         end_evalf

forino_evalf:

    vmovq   xmm2,[rax+rsi*8]
    vmovq   xmm3,xmm2
    vmulpd  xmm3,xmm3
    vaddpd  xmm0,xmm3

    vmovq   xmm3,[rbx+rsi*8]
    vmulpd  xmm2,xmm3
    vaddpd  xmm1,xmm2

    inc     rsi
    cmp     rsi,rdi
    jl      forino_evalf

end_evalf:
    vmovq [r8],xmm0
    vmovq [r9],xmm1

    stop
    



;compute_avg_64(MATRIX x, int np, int d, VECTOR c,type den, VECTOR ris)
section .data
section .bss
section .text

global compute_avg_64
    ; rdi = x 
    ; rsi = np
    ; rdx = d
    ; rcx = c
    ; r8 = den
    ; r9 = ris
    UNROLL_W equ 8
compute_avg_64:

    start
    mov rax,rdi     ; x
    mov rbx,rcx     ; c
    mov r10,rdx     ; d
    mov r11,rsi     ; np
    

    vbroadcastsd ymm0,[r8]
    
    mov rsi,0
    
    
foriw:
    
    mov             rdi,0
    vbroadcastsd    ymm7,[rbx+rsi*8]
    
    
    vdivpd  ymm7,ymm0
    
    mov     rdx,r10
    imul    rdx,rsi

    sub     r10,UNROLL_W-1
forjw:
    
    add     rdx,rdi
    vmovupd ymm1,[rax+rdx*8]
    sub     rdx,rdi
    vmulpd  ymm1,ymm7
    vaddpd  ymm1,[r9+rdi*8]
    vmovapd [r9+rdi*8],ymm1


    add     rdx,rdi
    vmovupd ymm1,[rax+rdx*8+32]
    sub     rdx,rdi
    vmulpd  ymm1,ymm7
    vaddpd  ymm1,[r9+rdi*8+32]
    vmovupd [r9+rdi*8+32],ymm1
    
    add     rdi,UNROLL_W
    cmp     rdi,r10
    jl      forjw

    add     r10,UNROLL_W-1

    cmp     rdi,r10
    jge     endforjw
    
forinow:
    add     rdx,rdi
    vmovq   xmm1,[rax+rdx*8]
    sub     rdx,rdi
    vmulsd  xmm1,xmm7
    vaddsd  xmm1,[r9+rdi*8]
    vmovq   [r9+rdi*8],xmm1
    
    inc rdi
    cmp rdi,r10
    jl forinow
    

endforjw:

    inc rsi 
    cmp rsi,r11
    jl  foriw


    stop
