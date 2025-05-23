# SPDX-License-Identifier: GPL-2.0
TCOBJ= tc.o tc_qdisc.o tc_class.o tc_filter.o tc_util.o tc_monitor.o \
       tc_exec.o m_police.o m_estimator.o m_action.o m_ematch.o \
       emp_ematch.tab.o emp_ematch.lex.o

include ../config.mk

SHARED_LIBS ?= y

TCMODULES :=
TCMODULES += q_fifo.o
TCMODULES += q_sfq.o
TCMODULES += q_red.o
TCMODULES += q_prio.o
TCMODULES += q_skbprio.o
TCMODULES += q_tbf.o
TCMODULES += q_multiq.o
TCMODULES += q_netem.o
TCMODULES += q_dlc.o
TCMODULES += q_choke.o
TCMODULES += q_sfb.o
TCMODULES += f_u32.o
TCMODULES += f_route.o
TCMODULES += f_fw.o
TCMODULES += f_basic.o
TCMODULES += f_bpf.o
TCMODULES += f_flow.o
TCMODULES += f_cgroup.o
TCMODULES += f_flower.o
TCMODULES += q_gred.o
TCMODULES += q_ingress.o
TCMODULES += q_hfsc.o
TCMODULES += q_htb.o
TCMODULES += q_drr.o
TCMODULES += q_qfq.o
TCMODULES += m_gact.o
TCMODULES += m_mirred.o
TCMODULES += m_mpls.o
TCMODULES += m_nat.o
TCMODULES += m_pedit.o
TCMODULES += m_ife.o
TCMODULES += m_skbedit.o
TCMODULES += m_skbmod.o
TCMODULES += m_csum.o
TCMODULES += m_simple.o
TCMODULES += m_vlan.o
TCMODULES += m_connmark.o
TCMODULES += m_ctinfo.o
TCMODULES += m_bpf.o
TCMODULES += m_tunnel_key.o
TCMODULES += m_sample.o
TCMODULES += m_ct.o
TCMODULES += m_gate.o
TCMODULES += p_ip.o
TCMODULES += p_ip6.o
TCMODULES += p_icmp.o
TCMODULES += p_eth.o
TCMODULES += p_tcp.o
TCMODULES += p_udp.o
TCMODULES += em_nbyte.o
TCMODULES += em_cmp.o
TCMODULES += em_u32.o
TCMODULES += em_canid.o
TCMODULES += em_meta.o
TCMODULES += q_mqprio.o
TCMODULES += q_codel.o
TCMODULES += q_fq_codel.o
TCMODULES += q_fq.o
TCMODULES += q_pie.o
TCMODULES += q_fq_pie.o
TCMODULES += q_cake.o
TCMODULES += q_hhf.o
TCMODULES += q_clsact.o
TCMODULES += e_bpf.o
TCMODULES += f_matchall.o
TCMODULES += q_cbs.o
TCMODULES += q_etf.o
TCMODULES += q_taprio.o
TCMODULES += q_plug.o
TCMODULES += q_ets.o

TCSO :=

ifneq ($(TC_CONFIG_NO_XT),y)
  ifeq ($(TC_CONFIG_XT),y)
    TCMODULES += em_ipt.o
    ifeq ($(TC_CONFIG_IPSET),y)
      TCMODULES += em_ipset.o
    endif
  endif
endif

TCOBJ += $(TCMODULES)
LDLIBS += -L. -lm

ifeq ($(SHARED_LIBS),y)
LDLIBS += -ldl
LDFLAGS += -Wl,-export-dynamic
endif

TCLIB := tc_core.o
TCLIB += tc_red.o
TCLIB += tc_estimator.o
TCLIB += tc_stab.o
TCLIB += tc_qevent.o

CFLAGS += -DCONFIG_GACT -DCONFIG_GACT_PROB
ifneq ($(IPT_LIB_DIR),)
	CFLAGS += -DIPT_LIB_DIR=\"$(IPT_LIB_DIR)\"
endif

LEX := flex
CFLAGS += -DYY_NO_INPUT

MODDESTDIR := $(DESTDIR)$(LIBDIR)/tc

%.so: %.c
	$(QUIET_CC)$(CC) $(CFLAGS) $(CPPFLAGS) $(LDFLAGS) -shared -fpic $< -o $@


all: tc $(TCSO)

tc: $(TCOBJ) $(LIBNETLINK) libtc.a
	$(QUIET_LINK)$(CC) $^ $(LDFLAGS) $(LDLIBS) -o $@

libtc.a: $(TCLIB)
	$(QUIET_AR)$(AR) rcs $@ $^

install: all
	mkdir -p $(MODDESTDIR)
	install -m 0755 tc $(DESTDIR)$(SBINDIR)
	for i in $(TCSO); \
	do install -m 755 $$i $(MODDESTDIR); \
	done
	if [ ! -f $(MODDESTDIR)/m_ipt.so ]; then \
	if [ -f $(MODDESTDIR)/m_xt.so ]; \
		then ln -s m_xt.so $(MODDESTDIR)/m_ipt.so ; \
	elif [ -f $(MODDESTDIR)/m_xt_old.so ]; \
		then ln -s m_xt_old.so $(MODDESTDIR)/m_ipt.so ; \
	fi; \
	fi

clean:
	rm -f $(TCOBJ) $(TCLIB) libtc.a tc *.so emp_ematch.tab.h; \
	rm -f emp_ematch.tab.*

m_xt.so: m_xt.c
	$(QUIET_CC)$(CC) $(CFLAGS) $(CPPFLAGS) $(LDFLAGS) -shared -fpic -o m_xt.so m_xt.c $$($(PKG_CONFIG) xtables --cflags --libs)

m_xt_old.so: m_xt_old.c
	$(QUIET_CC)$(CC) $(CFLAGS) $(CPPFLAGS) $(LDFLAGS) -shared -fpic -o m_xt_old.so m_xt_old.c $$($(PKG_CONFIG) xtables --cflags --libs)

em_ipset.o: CFLAGS += $$($(PKG_CONFIG) xtables --cflags)

em_ipt.o: CFLAGS += $$($(PKG_CONFIG) xtables --cflags)

ifeq ($(TC_CONFIG_XT),y)
  LDLIBS += $$($(PKG_CONFIG) xtables --libs)
endif

%.tab.c: %.y
	$(QUIET_YACC)$(YACC) $(YACCFLAGS) -p ematch_ -b $(basename $(basename $@)) $<

%.lex.c: %.l
	$(QUIET_LEX)$(LEX) $(LEXFLAGS) -o$@ $<

# our lexer includes the header from yacc, so make sure
# we don't attempt to compile it before the header has
# been generated as part of the yacc step.
emp_ematch.lex.o: emp_ematch.tab.c

ifneq ($(SHARED_LIBS),y)

tc: static-syms.o
static-syms.o: static-syms.h
static-syms.h: $(wildcard *.c)
	files="$^" ; \
	for s in `grep -B 3 '\<dlsym' $$files | sed -n '/snprintf/{s:.*"\([^"]*\)".*:\1:;s:%s::;p}'` ; do \
		sed -n '/'$$s'[^ ]* =/{s:.* \([^ ]*'$$s'[^ ]*\) .*:extern char \1[] __attribute__((weak)); if (!strcmp(sym, "\1")) return \1;:;p}' $$files ; \
	done > $@

endif
